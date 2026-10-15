import numpy as np
import pytest
from scipy.stats import multivariate_normal

import gratis


@pytest.mark.parametrize("covariance", [-np.eye(2), [[1, 2], [2, 1]],
                                        [[1, 0], [0, 0]], [[1, 1], [0, 1]],
                                        [[np.nan, 0], [0, 1]], [[np.inf, 0], [0, 1]]])
@pytest.mark.parametrize("sampling", [True, False])
def test_invalid_covariances_are_rejected(covariance, sampling):
    arguments = dict(means=np.zeros((2, 1)), sigmas=np.asarray(covariance)[:, :, None], weights=[1])
    with pytest.raises(ValueError, match="sigma matrices"):
        if sampling:
            gratis.rmixnorm(2, **arguments, rng=1)
        else:
            gratis.dmixnorm([0, 0], **arguments)


@pytest.mark.parametrize("variance", [0, -1, np.nan, np.inf])
def test_invalid_scalar_variances_are_rejected(variance):
    with pytest.raises(ValueError, match="sigma matrices"):
        gratis.rmixnorm(1, [0], [variance], [1], rng=1)
    with pytest.raises(ValueError, match="sigma matrices"):
        gratis.dmixnorm([0], [0], [variance], [1])


def test_zero_weight_does_not_hide_an_invalid_covariance():
    covariances = np.stack((np.eye(2), -np.eye(2)), axis=2)
    with pytest.raises(ValueError, match="positive definite"):
        gratis.rmixnorm(0, np.zeros((2, 2)), covariances, [1, 0], rng=1)
    with pytest.raises(ValueError, match="positive definite"):
        gratis.dmixnorm([0, 0], np.zeros((2, 2)), covariances, [1, 0])


def test_valid_density_and_sample_covariance_match_reference():
    mean = np.array([0.2, -0.1])
    covariance = np.array([[2, 0.3], [0.3, 1]])
    observations = np.array([[0, 0], [1, 2], [-2, 1]])
    actual = gratis.dmixnorm(observations, mean[:, None], covariance[:, :, None], [1], log=True)
    np.testing.assert_allclose(actual, multivariate_normal.logpdf(observations, mean, covariance))
    draws = gratis.rmixnorm(50000, mean[:, None], covariance[:, :, None], [1], rng=1)
    np.testing.assert_allclose(np.mean(draws, axis=0), mean, atol=0.03)
    np.testing.assert_allclose(np.cov(draws, rowvar=False), covariance, atol=0.05)


def test_log_density_handles_all_negative_infinite_components():
    assert gratis.dmixnorm([np.inf], [0], [1], [1], log=True)[0] == -np.inf
    assert gratis.dmixnorm([np.inf], [0], [1], [1])[0] == 0


@pytest.mark.parametrize("weights", [[np.nan, 1], [np.inf, 1], [0, 0], [-1, 1]])
def test_invalid_weights_are_rejected(weights):
    with pytest.raises(ValueError, match="weight"):
        gratis.rmixnorm(1, [0, 1], [1, 1], weights, rng=1)


def test_large_finite_weights_are_normalized_without_overflow():
    model = gratis.mar_model(k=2, p=0, d=0, constants=[0, 0], sigmas=[1, 1],
                             weights=[1e308, 1e308])
    np.testing.assert_allclose(model.weights, [0.5, 0.5])
