import importlib
from unittest.mock import patch

import numpy as np
import pytest
from scipy.signal import lfilter
from scipy.special import logsumexp
from scipy.stats import norm

import gratis


@pytest.mark.parametrize("arguments", [
    dict(order=(0, 0, 1), theta=[np.nan]),
    dict(order=(1, 0, 0), phi=[np.inf]),
    dict(phi=[.8]), dict(order=(1, 0, 0)), dict(order=(0, .5, 0)),
    dict(order=(0, 0)), dict(order=0), dict(order=None),
    dict(seasonal_order=(1, 0, 0), Phi=[.2]),
    dict(frequency=1.5), dict(frequency=0), dict(order=(0, 2, 0), constant=1),
    dict(order=(0, 2, 0), seasonal_order=(0, 1, 0), frequency=4),
    dict(sigma=np.nan), dict(sigma=0), dict(constant=np.inf),
])
def test_direct_arima_rejects_invalid_parameters(arguments):
    with pytest.raises(ValueError):
        gratis.ARIMAModel(**arguments)


@pytest.mark.parametrize("arguments", [
    dict(alpha=np.nan), dict(alpha=None), dict(alpha=0), dict(alpha=1),
    dict(trend="A", alpha=.2, beta=.3),
    dict(frequency=4, seasonal="A", alpha=.2, gamma=.9, season=[0]*4),
    dict(frequency=4, seasonal="M", season=[1, 1, 1, -1], level=1),
    dict(frequency=4, seasonal="A", season=[0]),
    dict(frequency=4, seasonal="A", season=[0, 0, np.nan, 0]),
    dict(seasonal="A"), dict(error="M", level=0), dict(level=np.nan),
    dict(slope=np.inf), dict(sigma=np.inf), dict(sigma=0), dict(frequency=1.5),
    dict(trend="A", damped=True, phi=0), dict(damped="yes"),
])
def test_direct_ets_rejects_invalid_parameters(arguments):
    with pytest.raises(ValueError):
        gratis.ETSModel(**arguments)


def test_direct_ets_checks_admissibility_in_addition_to_simple_bounds(monkeypatch):
    module = importlib.import_module("gratis.ets")
    monkeypatch.setattr(module, "ets_admissible", lambda *args: False)
    with pytest.raises(ValueError, match="inadmissible"):
        gratis.ETSModel(alpha=.2)


def test_direct_ets_accepts_numpy_scalar_smoothing_parameters():
    model = gratis.ETSModel(alpha=np.array(.2))
    np.testing.assert_allclose(model.simulate(5, rng=1), gratis.ETSModel(alpha=.2).simulate(5, rng=1))


@pytest.mark.parametrize("backend", ["numpy", "statsmodels"])
def test_direct_arima_agrees_with_factory_and_revalidates_mutated_coefficients(backend):
    parameters = dict(frequency=4, phi=[.3], theta=[.2], Phi=[.2], Theta=[-.1],
                      constant=1, sigma=.5, backend=backend)
    model = gratis.ARIMAModel(order=(1, 1, 1), seasonal_order=(1, 0, 1), **parameters)
    factory = gratis.arima_model(p=1, d=1, q=1, P=1, D=0, Q=1, **parameters)
    np.testing.assert_allclose(model.simulate(20, rng=1), factory.simulate(20, rng=1))
    model.theta[0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        model.simulate(5, rng=1)


def test_direct_ets_agrees_with_factory_and_copies_initial_states():
    states = np.array([-1., .5, 1., -.5])
    parameters = dict(frequency=4, error="A", trend="A", seasonal="A", damped=True,
                      alpha=.2, beta=.05, gamma=.1, phi=.9, level=10, slope=.2,
                      season=states, sigma=.5)
    model = gratis.ETSModel(**parameters)
    factory = gratis.ets_model(**parameters)
    np.testing.assert_allclose(model.simulate(20, rng=1), factory.simulate(20, rng=1))
    states[0] = 100
    np.testing.assert_array_equal(model.season, [-1, .5, 1, -.5])


def test_factory_and_constructor_share_cached_admissibility_work():
    module = importlib.import_module("gratis.ets")
    module.ets_admissible.cache_clear()
    with patch.object(np, "roots", wraps=np.roots) as roots:
        gratis.ets_model(frequency=4, error="A", trend="N", seasonal="A",
                         alpha=.2, gamma=.1, level=1, season=[0]*4, sigma=1)
    assert roots.call_count == 1


@pytest.mark.parametrize("name", ["P", "D", "Q"])
def test_nonseasonal_arima_does_not_truncate_fractional_seasonal_orders(name):
    with pytest.raises(ValueError):
        gratis.arima_model(p=0, d=0, q=0, **{name: .9})


def test_statsmodels_arima_uses_filter_without_dense_state_space(monkeypatch):
    from statsmodels.tsa import arima_process
    import statsmodels.tsa.statespace.sarimax as sarimax

    def unexpected_state_space(*args, **kwargs):
        raise AssertionError("dense seasonal state space must not be allocated")

    monkeypatch.setattr(sarimax, "SARIMAX", unexpected_state_space)
    model = gratis.arima_model(frequency=365, p=1, d=0, q=1, P=1, D=0, Q=1,
                               phi=[.3], theta=[.2], Phi=[.5], Theta=[.2],
                               constant=0, sigma=1)
    innovations = np.random.default_rng(1).normal(size=1500)
    expected = lfilter(model.ma_polynomial(), model.ar_polynomial(include_differences=False), innovations)
    with patch.object(arima_process, "arma_generate_sample", wraps=arima_process.arma_generate_sample) as sample:
        actual = model.simulate(innov=innovations, n_start=0)
    assert sample.call_count == 1
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_statsmodels_filter_does_not_silently_fall_back(monkeypatch):
    module = importlib.import_module("gratis.arima")
    monkeypatch.setattr(module, "_load_statsmodels_arma", lambda: (None, ImportError("unavailable")))
    with pytest.raises(ImportError, match="statsmodels is required"):
        gratis.ARIMAModel().simulate(5)


@pytest.mark.parametrize("n", [0, 1, 2, 3, 100])
def test_vectorized_conditional_density_matches_pointwise_reference(n):
    y = np.random.default_rng(1).normal(size=n)
    parameters = [[.1, .2], [-.1, .4, .1], [0, .6]]
    variances = np.array([np.linspace(.5, 2, n), np.ones(n), np.full(n, 2.)])
    weights = np.array([.2, .3, .5])
    expected = np.zeros(n)
    for t in range(2, n):
        means = [np.dot(p, np.r_[1., y[t-np.arange(1, len(p))]]) for p in parameters]
        expected[t] = logsumexp(np.log(weights) + norm.logpdf(y[t], means, np.sqrt(variances[:, t])))
    np.testing.assert_allclose(gratis.dmixnorm_ts(y, parameters, variances, weights, log=True), expected)
    np.testing.assert_allclose(gratis.dmixnorm_ts(y, parameters, variances, weights), np.exp(expected))


def test_inactive_components_are_not_evaluated_on_extreme_observations():
    with np.errstate(all="raise"):
        actual = gratis.dmixnorm_ts([1e308, 1e308, 0], [[0, 0, 0], [0, 2, -2]],
                                    [[1]*3, [1]*3], [1, 0], log=True)
    np.testing.assert_allclose(actual, [0, 0, norm.logpdf(0)])
