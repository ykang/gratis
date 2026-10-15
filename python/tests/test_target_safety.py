import importlib

import numpy as np
import pytest

import gratis

target_module = importlib.import_module("gratis.target")


def _parameters(weights):
    k = len(weights)
    return np.r_[np.zeros(6 * k), np.ones(k), np.cumsum(weights)[:-1]]


@pytest.mark.parametrize("weights", [[1], [0.9, 0.05, 0.05], [0.05, 0.9, 0.05],
                                    [0, 0, 1], [1, 0, 0], [0.2, 0, 0.8]])
def test_target_weights_cover_interior_and_boundary_of_simplex(weights):
    model = target_module._params_to_mar(_parameters(weights), np.array([1]), len(weights))
    np.testing.assert_allclose(model.weights, weights, atol=1e-15)
    assert np.isfinite(model.simulate(12, n_start=0, rng=1)).all()


@pytest.mark.parametrize("x", [[1e155, 2e155, 3e155], [1000001, 1000002, 1000003],
                               [-1e308, 0, 1e308], [1e-300, 2e-300, 3e-300],
                               [1, 1 + np.finfo(float).eps, 1 + 2 * np.finfo(float).eps]])
def test_scaling_is_stable_across_magnitude_and_offset(x):
    with np.errstate(all="raise"):
        scaled = gratis.scalets(x)
    np.testing.assert_allclose(scaled, [-1, 0, 1], atol=1e-12)
    np.testing.assert_allclose(np.std(scaled, ddof=1), 1)


@pytest.mark.parametrize("x", [[], [5], [1e155, 1e155], [0, 0, 0]])
def test_scaling_preserves_exact_constants_and_empty_inputs(x):
    np.testing.assert_array_equal(gratis.scalets(x), x)


@pytest.mark.parametrize("x", [[np.inf, 1], [np.nan, 1]])
def test_scaling_rejects_nonfinite_inputs(x):
    with pytest.raises(ValueError, match="finite"):
        gratis.scalets(x)


def test_explosive_finite_candidate_does_not_scale_to_zero():
    params = np.array([0, 0, 1.5, 0, 0, 0, 1.])
    score, series = target_module._fitness_mar(
        params, length=1000, seasonal_periods=np.array([1]), k=1,
        feature_function=lambda x: [np.std(x, ddof=1)], target=np.array([0.]),
        rng=np.random.default_rng(1),
    )
    assert np.isfinite(series).all()
    np.testing.assert_allclose(score, -1)
    np.testing.assert_allclose(np.std(series, ddof=1), 1)


def test_nonfinite_candidate_is_rejected_before_callback():
    def unexpected_callback(x):
        raise AssertionError("callback must not see nonfinite data")

    score, series = target_module._fitness_mar(
        np.array([0, 0, 1.5, 0, 0, 0, 1.]), length=4000,
        seasonal_periods=np.array([1]), k=1, feature_function=unexpected_callback,
        target=np.array([0.]), rng=np.random.default_rng(1),
    )
    assert score == -np.inf
    assert np.isnan(series).all()


def _search(**kwargs):
    arguments = dict(length=20, k=1, feature_function=lambda x: [0.], target=[0.],
                     population_size=4, max_iter=2, rng=1)
    arguments.update(kwargs)
    return gratis.simulate_target(**arguments)


def test_target_search_reports_unmet_tolerance_and_returns_best_series():
    with pytest.warns(RuntimeWarning, match="mean absolute feature error=1, tolerance=0.01"):
        result = _search(target=[1.], tolerance=0.01)
    assert result.shape == (20,)
    assert np.isfinite(result).all()


def test_target_search_supports_one_component_and_successful_tolerance():
    result = _search(feature_function=lambda x: [np.mean(x)], tolerance=1e-8)
    assert result.shape == (20,)
    assert abs(np.mean(result)) <= 1e-8


def test_target_callback_exceptions_are_not_suppressed():
    def broken(x):
        raise LookupError("feature callback bug")

    with pytest.raises(LookupError, match="feature callback bug"):
        _search(feature_function=broken)


def test_invalid_feature_lengths_and_all_invalid_features_are_reported():
    with pytest.raises(ValueError, match="one value per target feature"):
        _search(feature_function=lambda x: [0., 1.])
    with pytest.raises(RuntimeError, match="finite series"):
        _search(feature_function=lambda x: [np.nan])


@pytest.mark.parametrize("kwargs", [{"target": []}, {"target": [np.nan]},
                                     {"tolerance": -1}, {"tolerance": np.inf},
                                     {"seasonal_periods": 1.5}, {"length": 2.5}, {"k": 1.5}])
def test_invalid_target_search_inputs_raise(kwargs):
    with pytest.raises(ValueError):
        _search(**kwargs)


@pytest.mark.parametrize("n", [0, 1, 5])
def test_zero_order_mar_lengths_at_zero_burnin(n):
    model = gratis.mar_model(k=1, p=0, d=0, constants=0, sigmas=1, weights=[1])
    assert model.simulate(n, n_start=0, rng=1).shape == (n,)


@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_short_mixture_paths_and_densities_preserve_initialization(n):
    pars, sigmas = [[0, .5, .2, .1]], [np.ones(n)]
    series = gratis.rmixnorm_ts(n, pars, sigmas, [1], yinit=2, rng=1)
    np.testing.assert_array_equal(series, np.full(n, 2))
    np.testing.assert_array_equal(gratis.dmixnorm_ts(series, pars, sigmas, [1]), np.ones(n))
    np.testing.assert_array_equal(gratis.dmixnorm_ts(series, pars, sigmas, [1], log=True), np.zeros(n))
