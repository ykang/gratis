"""Mixture-normal simulation and density utilities."""

from __future__ import annotations

from collections.abc import Sequence
import numpy as np

from ._random import SeedLike, get_rng


def _normalise_weights(weights: np.ndarray, k: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or weights.size != k:
        raise ValueError("weights must be a one-dimensional vector of length k")
    if np.any(weights < 0):
        raise ValueError("weights must be non-negative")
    total = weights.sum()
    if total <= 0:
        raise ValueError("at least one mixture weight must be positive")
    return weights / total


def _coerce_mixture_inputs(
    means: np.ndarray | Sequence[float],
    sigmas: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray, bool]:
    means = np.asarray(means, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    univariate = means.ndim == 1

    if univariate:
        means = means.reshape(1, -1)
        sigmas = sigmas.reshape(1, 1, -1)
    elif means.ndim != 2:
        raise ValueError("means must be a vector or a q x k matrix")

    q, k = means.shape
    if sigmas.shape != (q, q, k):
        raise ValueError("sigmas must be a q x q x k array")

    return means, sigmas, univariate


def rmixnorm(
    n: int,
    means: np.ndarray | Sequence[float],
    sigmas: np.ndarray | Sequence[float],
    weights: Sequence[float],
    *,
    rng: SeedLike = None,
) -> np.ndarray:
    """Draw samples from a finite mixture of multivariate normal distributions.

    Parameters mirror the R function: `means` is a q x k matrix and `sigmas` is
    a q x q x k covariance array. For univariate mixtures, vectors are accepted.
    """
    if n < 0:
        raise ValueError("n must be non-negative")

    means, sigmas, univariate = _coerce_mixture_inputs(means, sigmas)
    q, k = means.shape
    weights = _normalise_weights(np.asarray(weights, dtype=float), k)
    generator = get_rng(rng)

    components = generator.choice(k, size=n, replace=True, p=weights)
    out = np.empty((n, q), dtype=float)
    for component in range(k):
        rows = components == component
        count = int(rows.sum())
        if count:
            out[rows, :] = generator.multivariate_normal(
                mean=means[:, component],
                cov=sigmas[:, :, component],
                size=count,
                check_valid="ignore",
            )

    return out[:, 0] if univariate else out


def _coerce_observations(x: np.ndarray | Sequence[float], q: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if q == 1:
        return x.reshape(-1, 1)
    if x.ndim == 1 and x.size == q:
        return x.reshape(1, q)
    if x.ndim == 2 and x.shape[1] == q:
        return x
    raise ValueError("x must be a vector of length q or an n x q matrix")


def _log_multivariate_normal_density(x: np.ndarray, mean: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    sign, logdet = np.linalg.slogdet(sigma)
    if sign <= 0:
        raise ValueError("sigma matrices must be positive definite")
    centered = x - mean
    solved = np.linalg.solve(sigma, centered.T).T
    quadratic = np.sum(centered * solved, axis=1)
    q = x.shape[1]
    return -0.5 * (q * np.log(2.0 * np.pi) + logdet + quadratic)


def _logsumexp(values: np.ndarray, axis: int) -> np.ndarray:
    maximum = np.max(values, axis=axis, keepdims=True)
    return (maximum + np.log(np.sum(np.exp(values - maximum), axis=axis, keepdims=True))).squeeze(axis)


def dmixnorm(
    x: np.ndarray | Sequence[float],
    means: np.ndarray | Sequence[float],
    sigmas: np.ndarray | Sequence[float],
    weights: Sequence[float],
    *,
    log: bool = False,
) -> np.ndarray:
    """Evaluate the density of a finite mixture of multivariate normals."""
    means, sigmas, univariate = _coerce_mixture_inputs(means, sigmas)
    q, k = means.shape
    weights = _normalise_weights(np.asarray(weights, dtype=float), k)
    observations = _coerce_observations(x, q)

    component_log_density = np.empty((observations.shape[0], k), dtype=float)
    for component in range(k):
        if weights[component] == 0:
            component_log_density[:, component] = -np.inf
        else:
            component_log_density[:, component] = (
                np.log(weights[component])
                + _log_multivariate_normal_density(observations, means[:, component], sigmas[:, :, component])
            )
    log_density = _logsumexp(component_log_density, axis=1)
    density = log_density if log else np.exp(log_density)
    if univariate and np.asarray(x).ndim == 0:
        return np.asarray(density[0])
    return density


def rmixnorm_ts(
    n: int,
    means_ar_par_list: Sequence[Sequence[float]],
    sigmas_list: Sequence[Sequence[float]],
    weights: Sequence[float],
    *,
    yinit: float = 0.0,
    rng: SeedLike = None,
) -> np.ndarray:
    """Simulate an autoregressive process with mixture-normal innovations."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if len(means_ar_par_list) == 0:
        raise ValueError("means_ar_par_list must contain at least one component")
    if len(sigmas_list) != len(means_ar_par_list):
        raise ValueError("sigmas_list must match means_ar_par_list")

    ar_params = [np.asarray(params, dtype=float) for params in means_ar_par_list]
    n_lags = np.asarray([params.size - 1 for params in ar_params], dtype=int)
    if np.any(n_lags < 1):
        raise ValueError("each AR parameter vector must include an intercept and at least one lag")

    k = len(ar_params)
    weights = _normalise_weights(np.asarray(weights, dtype=float), k)
    sigmas = np.column_stack([np.asarray(values, dtype=float) for values in sigmas_list])
    if sigmas.shape != (n, k):
        raise ValueError("each sigma vector must have length n")

    generator = get_rng(rng)
    y = np.full(n, yinit, dtype=float)
    max_lag = int(n_lags.max())
    for index in range(max_lag, n):
        means = np.empty(k, dtype=float)
        for component, params in enumerate(ar_params):
            lag_count = n_lags[component]
            previous = y[index - np.arange(1, lag_count + 1)]
            means[component] = np.dot(params, np.concatenate(([1.0], previous)))
        y[index] = rmixnorm(
            1,
            means=means,
            sigmas=sigmas[index, :],
            weights=weights,
            rng=generator,
        )[0]

    return y


def dmixnorm_ts(
    y: Sequence[float] | np.ndarray,
    means_ar_par_list: Sequence[Sequence[float]],
    sigmas_list: Sequence[Sequence[float]],
    weights: Sequence[float],
    *,
    log: bool = False,
) -> np.ndarray:
    """Evaluate conditional mixture-normal densities for an AR mixture series."""
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.size
    if len(means_ar_par_list) == 0:
        raise ValueError("means_ar_par_list must contain at least one component")
    if len(sigmas_list) != len(means_ar_par_list):
        raise ValueError("sigmas_list must match means_ar_par_list")

    ar_params = [np.asarray(params, dtype=float) for params in means_ar_par_list]
    n_lags = np.asarray([params.size - 1 for params in ar_params], dtype=int)
    if np.any(n_lags < 1):
        raise ValueError("each AR parameter vector must include an intercept and at least one lag")

    k = len(ar_params)
    weights = _normalise_weights(np.asarray(weights, dtype=float), k)
    sigmas = np.column_stack([np.asarray(values, dtype=float) for values in sigmas_list])
    if sigmas.shape != (n, k):
        raise ValueError("each sigma vector must have length len(y)")

    max_lag = int(n_lags.max())
    out_log = np.zeros(n, dtype=float)
    for index in range(max_lag, n):
        means = np.empty(k, dtype=float)
        for component, params in enumerate(ar_params):
            lag_count = n_lags[component]
            previous = y[index - np.arange(1, lag_count + 1)]
            means[component] = np.dot(params, np.concatenate(([1.0], previous)))
        out_log[index] = float(
            dmixnorm(
                y[index],
                means=means,
                sigmas=sigmas[index, :],
                weights=weights,
                log=True,
            )
        )

    return out_log if log else np.exp(out_log)
