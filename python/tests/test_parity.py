import numpy as np
import pytest

import gratis


def test_pi_coefficients_matches_r_for_seasonal_ar():
    expected = np.array([0.30, 0.00, 0.00, 0.20, -0.06])

    result = gratis.pi_coefficients(ar=[0.3], sar=[0.2], d=0, D=0, m=4)

    np.testing.assert_allclose(result, expected)


def test_pi_coefficients_matches_r_for_differenced_seasonal_ar():
    expected = np.array([1.30, -0.30, 0.00, 1.20, -1.56, 0.36, 0.00, -0.20, 0.26, -0.06])

    result = gratis.pi_coefficients(ar=[0.3], sar=[0.2], d=1, D=1, m=4)

    np.testing.assert_allclose(result, expected)


def test_pi_coefficients_matches_r_for_arima_with_ma():
    expected = np.array([1.30, -0.47, 0.153, 0.0153, 0.00153, 0.000153, 0.0000153, 0.00000153, 0.000000153])

    result = gratis.pi_coefficients(ar=[0.4, -0.2], ma=[0.1], d=1)

    np.testing.assert_allclose(result, expected)


def test_mar_model_matches_r_structure_for_fixed_inputs():
    model = gratis.mar_model(
        phi=np.array([[0.8, 0.6], [0.0, 0.3]]),
        d=0,
        constants=[1, -1],
        sigmas=[1, 2],
        weights=[0.8, 0.2],
    )

    np.testing.assert_array_equal(model.p, np.array([2, 2]))
    np.testing.assert_array_equal(model.d, np.array([0, 0]))
    np.testing.assert_allclose(model.constants, np.array([1, -1]))
    np.testing.assert_allclose(model.weights, np.array([0.8, 0.2]))
    np.testing.assert_allclose(model.ar, np.array([[1.0, -1.0], [0.8, 0.6], [0.0, 0.3]]))


def test_mar_model_matches_r_zero_order_zero_constant_structure():
    model = gratis.mar_model(k=1, p=0, d=0, constants=0, sigmas=1, weights=1)

    np.testing.assert_allclose(model.ar, np.array([[0.0]]))


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"p": 2, "phi": [0.1]}, "dimension of phi"),
        ({"frequency": 1, "P": 1}, "seasonal arguments"),
        ({"p": -1}, "p must be non-negative"),
    ],
)
def test_arima_model_validation_matches_r_style(kwargs, message):
    with pytest.raises(ValueError, match=message):
        gratis.arima_model(**kwargs)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"k": 2, "weights": [1]}, "dimension of weights"),
        ({"k": 2, "weights": [1, 0]}, "weights must be positive"),
        ({"k": 2, "seasonal_periods": [12, 24, 168]}, "dimension of seasonal_periods"),
    ],
)
def test_mar_model_validation_matches_r_style(kwargs, message):
    with pytest.raises(ValueError, match=message):
        gratis.mar_model(**kwargs)
