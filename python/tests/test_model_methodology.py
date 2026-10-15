import importlib

import numpy as np
import pytest
from statsmodels.tsa.arima_process import arma2ar
from scipy.signal import lfilter

import gratis


def test_ma_expansion_matches_independent_statsmodels_reference():
    expected = -arma2ar([1, -0.3], [1, 0.5], lags=501)[1:]
    expected = expected[:np.flatnonzero(abs(expected) > 1e-7)[-1] + 1]
    np.testing.assert_allclose(gratis.pi_coefficients(ar=0.3, ma=0.5), expected)
    assert gratis.pi_coefficients(ar=0.8, ma=-0.8).size == 0
    with pytest.warns(RuntimeWarning, match="increase n_terms"):
        gratis.pi_coefficients(ma=0.99, n_terms=10)


def test_annual_mar_preserves_all_seasonal_lags():
    expected = np.zeros(731)
    expected[np.array([1, 365, 366, 730, 731]) - 1] = [0.3, 0.2, -0.06, 0.3, -0.09]
    np.testing.assert_allclose(gratis.pi_coefficients(ar=0.3, sar=[0.2, 0.3], m=365), expected)
    model = gratis.mar_model(phi=[0.3], Phi=[0.2, 0.3], d=0, D=0, constants=[0],
                             sigmas=[1], weights=[1], seasonal_periods=365)
    np.testing.assert_allclose(model.ar[:, 0], np.r_[0, expected])


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
@pytest.mark.parametrize("d,D", [(0, 0), (1, 0), (0, 1)])
def test_arima_constant_is_mean_or_linear_drift(backend, d, D):
    model = gratis.arima_model(frequency=4, p=1, d=d, q=0, P=1, D=D, Q=0,
                               phi=[0.5], Phi=[0.2], constant=2, sigma=1e-10,
                               backend=backend)
    expected = np.full(12, 2.) if d + D == 0 else 2 * np.arange(1, 13)
    np.testing.assert_allclose(model.simulate(12, rng=1), expected, atol=1e-7)
    with pytest.raises(ValueError, match="must be zero"):
        gratis.arima_model(p=0, d=2, q=0, constant=1)


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
def test_arima_scales_burnin_and_output_together(backend):
    parameters = dict(p=1, d=0, q=1, phi=[0.99], theta=[0.3], constant=0, backend=backend)
    unit = gratis.arima_model(**parameters, sigma=1)
    tiny = gratis.arima_model(**parameters, sigma=1e-8)
    x = unit.simulate(100, rng=1)
    y = tiny.simulate(100, rng=1)
    np.testing.assert_allclose(y, x * 1e-8, rtol=1e-12, atol=1e-20)
    assert np.max(np.abs(y)) < 1e-5


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
def test_arima_automatic_burnin_uses_near_unit_root_decay(backend):
    phi, sigma = 0.999, 0.01
    model = gratis.arima_model(p=1, d=0, q=0, phi=[phi], constant=0, sigma=sigma,
                               backend=backend)
    burn = model._default_burnin()
    assert phi ** (2 * (burn + 1)) < 1e-8
    noise = np.random.default_rng(1).normal(0, sigma, burn + 3)
    expected = np.array([np.dot(noise[:t], phi ** np.arange(t - 1, -1, -1))
                         for t in range(burn + 1, burn + 4)])
    np.testing.assert_allclose(model.simulate(3, rng=1), expected, rtol=1e-12, atol=1e-12)
    assert model._default_burnin() == burn


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
@pytest.mark.parametrize("component", ["AR", "MA"])
def test_arima_preserves_zero_sum_seasonal_coefficients(backend, component):
    model = gratis.arima_model(frequency=4, p=0, d=0, q=0, P=2 if component == "AR" else 0,
                               D=0, Q=2 if component == "MA" else 0,
                               Phi=[0.5, -0.5] if component == "AR" else [],
                               Theta=[0.5, -0.5] if component == "MA" else [],
                               constant=0, sigma=1, backend=backend)
    innovations = np.r_[1., np.zeros(11)]
    expected = np.zeros(12)
    expected[[0, 4, 8]] = [1, 0.5, -0.25 if component == "AR" else -0.5]
    np.testing.assert_allclose(model.simulate(innov=innovations, n_start=0), expected)
    np.testing.assert_allclose(model.simulate(innov=[1.], n_start=0), [1.])


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
def test_arima_controlled_innovations_match_independent_filter(backend):
    model = gratis.arima_model(p=1, d=0, q=1, phi=[0.5], theta=[0.2], constant=0,
                               sigma=2, backend=backend)
    innovations = np.array([0.1, -0.3, 0.2, 0.5])
    expected = lfilter([1, 0.2], [1, -0.5], innovations)
    np.testing.assert_allclose(model.simulate(innov=innovations, n_start=0), expected)
    expected = lfilter([1, 0.2], [1, -0.5], np.r_[0.4, innovations])[1:]
    np.testing.assert_allclose(model.simulate(innov=innovations, start_innov=[0.4]), expected)
    for kwargs, message in [({"n_start": -1}, "non-negative"),
                            ({"n_start": 1.5}, "integer"),
                            ({"innov": [np.nan]}, "finite numeric vector"),
                            ({"start_innov": [np.inf]}, "finite numeric vector"),
                            ({"n_start": 2, "start_innov": [1]}, "length n_start")]:
        with pytest.raises(ValueError, match=message):
            model.simulate(**kwargs)


@pytest.mark.parametrize("kwargs", [{"alpha": 3}, {"alpha": 0.2, "beta": 0.3},
                                    {"phi": 0, "damped": True}, {"alpha": np.nan}])
def test_ets_rejects_invalid_fixed_parameters(kwargs):
    with pytest.raises(ValueError):
        gratis.ets_model(error="A", trend="A", seasonal="N", **kwargs)
    assert not gratis.ets_admissible(0.2, 0.1, phi=0)


def test_ets_resamples_after_rejection_and_bounds_retries(monkeypatch):
    module = importlib.import_module("gratis.ets")
    draws = []

    def check(alpha, *args):
        draws.append(alpha)
        return len(draws) > 1

    monkeypatch.setattr(module, "ets_admissible", check)
    model = gratis.ets_model(error="A", trend="N", seasonal="N", rng=1)
    assert len(draws) == 3 and draws[0] != draws[1]
    assert draws[2] == draws[1]  # The public model constructor validates the accepted draw.
    assert model.alpha == draws[1]
    monkeypatch.setattr(module, "ets_admissible", lambda *args: False)
    with pytest.raises(ValueError, match="inadmissible"):
        gratis.ets_model(error="A", trend="N", seasonal="N", alpha=0.2)
    with pytest.raises(ValueError, match="could not find admissible"):
        gratis.ets_model(error="A", trend="N", seasonal="N")


@pytest.mark.parametrize("period", [4, 12, 365])
def test_multiplicative_seasonal_factors_are_positive_and_normalized(period):
    model = gratis.ets_model(frequency=period, error="M", trend="N", seasonal="M",
                             alpha=0.2, gamma=0.1, rng=0)
    assert np.all(model.season > 0)
    np.testing.assert_allclose(model.season.mean(), 1)
    with pytest.raises(ValueError, match="factors must be positive"):
        gratis.ets_model(frequency=4, seasonal="M", season=[1, 1, 1, -1])
    with pytest.raises(ValueError, match="positive initial level"):
        gratis.ets_model(error="M", level=0)


def test_ets_has_one_backend_and_does_not_silently_fallback(monkeypatch):
    model = gratis.ets_model(error="A", trend="N", seasonal="N", alpha=0.2)
    assert model.backend == "statsmodels"
    with pytest.raises(ValueError, match="ETS uses statsmodels"):
        gratis.ets_model(backend="numpy")
    import sys
    monkeypatch.setitem(sys.modules, "statsmodels.tsa.exponential_smoothing.ets", None)
    with pytest.raises(ImportError):
        model.simulate(5)
