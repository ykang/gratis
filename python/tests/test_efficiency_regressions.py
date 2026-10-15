import importlib
import threading
from unittest.mock import patch

import numpy as np
import pytest

import gratis


@pytest.mark.parametrize("argument", ["phi", "Phi", "constants", "sigmas"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_mar_rejects_nonfinite_parameters(argument, value):
    arguments = dict(k=1, p=0, d=0, constants=[0], sigmas=[1], weights=[1])
    arguments[argument] = [value]
    with pytest.raises(ValueError, match="finite"):
        gratis.mar_model(**arguments)


@pytest.mark.parametrize("argument", ["p", "d", "P", "D", "k", "seasonal_periods"])
@pytest.mark.parametrize("value", [1.9, np.nan, np.inf, -1])
def test_mar_rejects_invalid_orders_instead_of_truncating(argument, value):
    arguments = dict(k=1, p=0, d=0, constants=[0], sigmas=[1], weights=[1])
    arguments[argument] = value
    with pytest.raises(ValueError):
        gratis.mar_model(**arguments)


@pytest.mark.parametrize("argument", ["ar", "ma", "sar", "sma"])
def test_polynomial_expansion_rejects_nan(argument):
    with pytest.raises(ValueError, match="finite"):
        gratis.pi_coefficients(**{argument: [np.nan]})


def dense_mar_reference(model, n, seed):
    rng = np.random.default_rng(seed)
    order = model.ar.shape[0] - 1
    total = n + 100 + order
    components = rng.choice(model.k, size=total, p=model.weights)
    noise = rng.normal(0, np.sqrt(model.sigmas[components]), total)
    y = np.zeros(total)
    for t in range(order, total):
        coefficients = model.ar[:, components[t]]
        y[t] = coefficients[0] + np.dot(coefficients[1:], y[t - np.arange(1, order + 1)]) + noise[t]
    return y[-n:]


@pytest.mark.parametrize("period", [1, 12, 365])
def test_sparse_mar_matches_dense_recursion(period):
    model = gratis.mar_model(k=2, p=1, phi=[[.2, .1]], d=0, P=2, D=0,
                             Phi=[[.3, .5], [.1, -.1]], constants=[0, 0],
                             sigmas=[1, 2], weights=[.4, .6], seasonal_periods=period)
    np.testing.assert_allclose(model.simulate(300, rng=7), dense_mar_reference(model, 300, 7),
                               rtol=1e-12, atol=1e-12)


def mixture_ts_reference(n, parameters, variances, weights, seed):
    rng = np.random.default_rng(seed)
    y = np.zeros(n)
    for t in range(max(len(p) - 1 for p in parameters), n):
        means = [np.dot(p, np.r_[1., y[t - np.arange(1, len(p))]]) for p in parameters]
        y[t] = gratis.rmixnorm(1, means, np.asarray(variances)[:, t], weights, rng=rng)[0]
    return y


def test_scalar_mixture_fast_path_matches_general_sampler_and_density():
    parameters = [[.1, .2], [-.1, .4, .1], [0, .6]]
    variances = [np.linspace(.5, 2, 100), np.ones(100), np.full(100, 2)]
    weights = [.2, .3, .5]
    y = gratis.rmixnorm_ts(100, parameters, variances, weights, rng=7)
    np.testing.assert_allclose(y, mixture_ts_reference(100, parameters, variances, weights, 7),
                               rtol=1e-12, atol=1e-12)
    expected = np.zeros(100)
    for t in range(2, 100):
        means = [np.dot(p, np.r_[1., y[t - np.arange(1, len(p))]]) for p in parameters]
        expected[t] = gratis.dmixnorm(y[t], means, np.asarray(variances)[:, t], weights, log=True)
    np.testing.assert_allclose(gratis.dmixnorm_ts(y, parameters, variances, weights, log=True), expected)


@pytest.mark.parametrize("variance", [0, -1, np.nan, np.inf])
def test_scalar_fast_paths_validate_all_variances_even_without_recursion(variance):
    with pytest.raises(ValueError, match="positive and finite"):
        gratis.rmixnorm_ts(1, [[0, .2], [0, .5]], [[1], [variance]], [1, 0])
    with pytest.raises(ValueError, match="positive and finite"):
        gratis.dmixnorm_ts([0], [[0, .2], [0, .5]], [[1], [variance]], [1, 0])


def test_target_retains_valid_incumbent_after_invalid_generation(monkeypatch):
    module = importlib.import_module("gratis.target")
    calls = []

    def fitness(params, **kwargs):
        calls.append(1)
        return (-1., np.arange(20.)) if len(calls) <= 4 else (-np.inf, np.full(20, np.nan))

    monkeypatch.setattr(module, "_fitness_mar", fitness)
    with pytest.warns(RuntimeWarning, match="error=1"):
        result = gratis.simulate_target(length=20, k=1, feature_function=lambda x: [0.],
                                        target=[1.], tolerance=.01, population_size=4,
                                        max_iter=4, rng=3)
    assert len(calls) == 8
    np.testing.assert_array_equal(result, np.arange(20.))


def test_parallel_target_uses_threads_and_matches_serial_result():
    threads = set()
    lock = threading.Lock()

    def feature(x):
        with lock:
            threads.add(threading.get_ident())
        return [0.]

    arguments = dict(length=20, k=1, feature_function=feature, target=[1.],
                     population_size=8, max_iter=2, rng=3, tolerance=.01)
    with pytest.warns(RuntimeWarning):
        serial = gratis.simulate_target(**arguments)
    assert threads == {threading.get_ident()}
    threads.clear()
    with pytest.warns(RuntimeWarning):
        parallel = gratis.simulate_target(**arguments, parallel=True)
    assert threads and threading.get_ident() not in threads
    np.testing.assert_array_equal(serial, parallel)


def test_parallel_target_propagates_callback_exception():
    def broken(x):
        raise LookupError("callback failed")

    with pytest.raises(LookupError, match="callback failed"):
        gratis.simulate_target(length=20, k=1, feature_function=broken, target=[0.],
                               population_size=4, max_iter=1, rng=1, parallel=True)


@pytest.mark.parametrize("kind", ["arima", "ets"])
def test_batch_prepares_statsmodels_once_and_preserves_seeded_paths(kind):
    model = (gratis.arima_model(frequency=4, p=1, d=1, q=1, P=1, D=0, Q=1,
                                phi=[.3], theta=[.2], Phi=[.2], Theta=[.1], sigma=.5)
             if kind == "arima" else
             gratis.ets_model(frequency=4, trend="A", seasonal="A", damped=True,
                              alpha=.2, beta=.05, gamma=.1, phi=.9, level=10, slope=.2,
                              season=[-1, .5, 1, -.5], sigma=.5))
    rng = np.random.default_rng(7)
    expected = np.column_stack([model.simulate(30, rng=rng) for _ in range(4)])
    with patch.object(type(model), "_prepare_statsmodels", autospec=True,
                      side_effect=type(model)._prepare_statsmodels) as prepare:
        actual = gratis.generate(model, length=30, nseries=4, rng=7)
    assert prepare.call_count == 1
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_batch_preparation_does_not_cache_mutable_coefficients_between_calls():
    model = gratis.arima_model(p=1, d=0, q=0, phi=[.2], sigma=1)
    before = gratis.generate(model, length=20, nseries=2, rng=7)
    model.phi[0] = .8
    after = gratis.generate(model, length=20, nseries=2, rng=7)
    assert not np.array_equal(before, after)
    rng = np.random.default_rng(7)
    np.testing.assert_allclose(after, np.column_stack([model.simulate(20, rng=rng) for _ in range(2)]))


def test_batch_arima_preserves_explicit_innovations_and_zero_length():
    model = gratis.arima_model(p=1, d=0, q=0, phi=[.3], sigma=1)
    innovations = np.linspace(-1, 1, 20)
    expected = model.simulate(innov=innovations, n_start=0)
    actual = gratis.generate(model, length=20, nseries=3, innov=innovations, n_start=0)
    np.testing.assert_allclose(actual, np.tile(expected[:, None], (1, 3)))
    assert gratis.generate(model, length=0, nseries=3).shape == (0, 3)
