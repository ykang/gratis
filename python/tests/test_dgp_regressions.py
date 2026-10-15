"""Regression checks for initialization, parameter ownership and DGP moments."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.signal import lfilter

import gratis
from gratis.arima import stationary_ar


@pytest.mark.parametrize("seasonal", [False, True])
def test_repeated_pole_burnin_captures_stationary_variance(seasonal):
    phi = -np.polynomial.polynomial.polypow([1, -.9], 6)[1:]
    args = (dict(frequency=4, p=0, P=6, Phi=phi) if seasonal
            else dict(p=6, phi=phi))
    model = gratis.arima_model(**args, d=0, q=0, D=0, Q=0, constant=0, sigma=1)
    burn = model._default_burnin()
    impulse = lfilter([1], model.ar_polynomial(include_differences=False),
                      np.r_[1., np.zeros(8191)])
    ratio = np.sum(impulse[:burn + 1] ** 2) / np.sum(impulse ** 2)
    assert ratio > 1 - 1e-8
    assert burn >= 128


def test_explicit_burnin_skips_automatic_search(monkeypatch):
    model = gratis.arima_model(p=1, d=0, q=0, phi=[.5], constant=0, sigma=1)
    monkeypatch.setattr(type(model), "_default_burnin", lambda self: pytest.fail("unexpected search"))
    np.testing.assert_allclose(model.simulate(innov=[1., 0.], n_start=0), [1., .5])
    np.testing.assert_allclose(model.simulate(innov=[1.], start_innov=[2.]), [2.])
    expected = model.simulate(5, n_start=0, rng=7)
    actual = gratis.generate(model, length=5, nseries=1, n_start=0, rng=7)
    np.testing.assert_allclose(actual[:, 0], expected)
    model.phi[0] = 1.1
    with pytest.raises(ValueError, match="stationary"):
        model.simulate(2, n_start=0)


def test_batch_prepares_automatic_burnin_only_once(monkeypatch):
    model = gratis.arima_model(p=1, d=0, q=0, phi=[.5], constant=0, sigma=1)
    original = type(model)._default_burnin
    calls = []

    def burnin(self):
        calls.append(1)
        return original(self)

    monkeypatch.setattr(type(model), "_default_burnin", burnin)
    gratis.generate(model, length=5, nseries=3, rng=7)
    assert len(calls) == 1


def test_unconverged_automatic_burnin_fails_explicitly():
    model = gratis.arima_model(p=1, d=0, q=0, phi=[.99999999], constant=0, sigma=1)
    with pytest.raises(ValueError, match="supply n_start explicitly"):
        model._default_burnin()
    assert np.all(np.isfinite(model.simulate(2, n_start=10, rng=1)))


def test_nearly_cancelled_poles_do_not_hide_slow_impulse_energy():
    phi, residual = .9999, 1e-5
    model = gratis.arima_model(p=1, d=0, q=1, phi=[phi], theta=[residual - phi],
                               constant=0, sigma=1)
    burn = model._default_burnin()
    variance = 1 + residual**2 / (1 - phi**2)
    retained = 1 + residual**2 * (1 - phi**(2*burn)) / (1 - phi**2)
    assert retained / variance > 1 - 1e-8


@pytest.mark.parametrize("p", [3, 4, 5, 6])
def test_stationary_sampler_covers_every_real_root_count(p):
    counts = set()
    for seed in range(100):
        phi = stationary_ar(p, rng=seed)
        roots = np.roots(np.r_[1., -phi])
        assert np.all(np.abs(roots) < 1)
        counts.add(np.count_nonzero(np.abs(roots.imag) < 1e-7))
    assert counts == set(range(p % 2, p + 1, 2))


def test_mar_owns_readonly_parameter_arrays():
    parameters = dict(phi=np.array([[.2]]), Phi=np.array([[.3]]),
                      constants=np.array([.1]), sigmas=np.array([1.]),
                      weights=np.array([1.]), seasonal_periods=np.array([4]))
    model = gratis.mar_model(**parameters, d=0, D=0)
    expected = model.simulate(20, rng=7)
    for values in parameters.values():
        values[...] = 0
    assert model.phi[0, 0] == .2
    assert model.Phi[0, 0] == .3
    np.testing.assert_array_equal(model.simulate(20, rng=7), expected)
    for name in ("p", "d", "phi", "P", "D", "Phi", "constants", "ar", "sigmas", "weights", "m"):
        assert not getattr(model, name).flags.writeable
        with pytest.raises(ValueError, match="read-only"):
            getattr(model, name).flat[0] = 0


@pytest.mark.parametrize("sigma", [1e155, 1e-200])
@pytest.mark.parametrize("kind", ["arima", "ets"])
def test_noise_variance_range_is_rejected_consistently(kind, sigma):
    factory = gratis.arima_model if kind == "arima" else gratis.ets_model
    with pytest.raises(ValueError, match="normal variance"):
        factory(sigma=sigma)
    model = (factory(p=0, d=0, q=0, constant=0, sigma=1) if kind == "arima" else
             factory(error="A", trend="N", seasonal="N", alpha=.2, level=0, sigma=1))
    with pytest.raises(ValueError, match="normal variance"):
        replace(model, sigma=sigma)


@pytest.mark.parametrize("sigma", [1e150, 1e-150])
@pytest.mark.parametrize("kind", ["arima", "ets"])
def test_large_and_small_supported_scales_preserve_noise(kind, sigma):
    if kind == "arima":
        arguments = dict(p=0, d=0, q=0, constant=0)
        factory = gratis.arima_model
    else:
        arguments = dict(error="A", trend="N", seasonal="N", alpha=.2, level=0)
        factory = gratis.ets_model
    unit = factory(**arguments, sigma=1)
    model = factory(**arguments, sigma=sigma)
    np.testing.assert_allclose(model.simulate(20, rng=7) / sigma, unit.simulate(20, rng=7),
                               rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("order,k", [(0, 1), (0, 3), (1, 1), (8, 1)])
@pytest.mark.parametrize("n,n_start", [(1, 0), (5, 0), (30, 10)])
def test_mar_fast_paths_preserve_recursion_and_rng(order, k, n, n_start):
    coefficients = np.zeros((order, k))
    if order:
        coefficients[-1] = .4
    model = gratis.mar_model(k=k, phi=coefficients, d=0, constants=np.arange(k) + .2,
                             sigmas=np.arange(k) + 1., weights=np.ones(k))
    rng = np.random.default_rng(7)
    total = n + n_start + order
    components = rng.choice(k, size=total, p=model.weights)
    noise = rng.normal(0, np.sqrt(model.sigmas[components]), total)
    expected = np.zeros(total)
    for t in range(order, total):
        j = components[t]
        expected[t] = (model.ar[0, j] +
                       np.dot(model.ar[1:, j], expected[t - np.arange(1, order + 1)]) + noise[t])
    actual_rng = np.random.default_rng(7)
    np.testing.assert_allclose(model.simulate(n, n_start=n_start, rng=actual_rng), expected[-n:],
                               rtol=1e-12, atol=1e-12)
    assert actual_rng.random() == rng.random()


@pytest.mark.parametrize("mixture", [False, True])
def test_mar_monte_carlo_mean_variance_and_lag_one(mixture):
    a = np.array([.4, -.2]) if mixture else np.array([.6])
    c = np.array([.5, -.2]) if mixture else np.array([.2])
    v = np.array([1., 2.]) if mixture else np.array([1.])
    w = np.array([.7, .3]) if mixture else np.array([1.])
    mean = np.dot(w, c) / (1 - np.dot(w, a))
    second = (np.dot(w, c*c + v) + 2*np.dot(w, c*a)*mean) / (1 - np.dot(w, a*a))
    expected = np.array([mean, second - mean**2, np.dot(w, a)])
    model = gratis.mar_model(phi=a[None, :], d=0, constants=c, sigmas=v, weights=w)
    estimates = []
    for seed in range(20):
        x = model.simulate(2000, n_start=300, rng=seed)
        estimates.append([x.mean(), x.var(ddof=1), np.corrcoef(x[:-1], x[1:])[0, 1]])
    estimates = np.array(estimates)
    mcse = estimates.std(axis=0, ddof=1) / np.sqrt(len(estimates))
    assert np.all(np.abs(estimates.mean(axis=0) - expected) < 6*mcse + .005)
