"""Mixture-normal simulation and density utilities."""

from __future__ import annotations

from collections.abc import Sequence
import numpy as np

from ._random import SeedLike, get_rng
from .arima import _validate_order


def _normalise_weights(weights: np.ndarray, k: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or weights.size != k:
        raise ValueError("weights must be a one-dimensional vector of length k")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative")
    if not np.any(weights > 0):
        raise ValueError("at least one mixture weight must be positive")
    weights = weights / np.max(weights)
    return weights / weights.sum()


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
    if q == 0 or k == 0 or np.any(~np.isfinite(means)):
        raise ValueError("means must be nonempty and finite")
    covariances = np.moveaxis(sigmas, 2, 0)
    if np.any(~np.isfinite(covariances)) or not np.allclose(
        covariances, covariances.swapaxes(1, 2), rtol=1e-12, atol=0.0
    ):
        raise ValueError("sigma matrices must be finite, symmetric and positive definite")
    try:
        factors = np.linalg.cholesky(covariances)
    except np.linalg.LinAlgError as exc:
        raise ValueError("sigma matrices must be positive definite") from exc

    return means, factors, univariate


def rmixnorm(
    n: int,
    means: np.ndarray | Sequence[float],
    sigmas: np.ndarray | Sequence[float],
    weights: Sequence[float],
    *,
    rng: SeedLike = None,
) -> np.ndarray:
    """Draw samples from a finite mixture of multivariate normal distributions.

    Parameters
    ----------
    n : int
        Non-negative number of observations.
    means : array_like, shape (k,) or (q, k)
        Finite component means. Columns identify components; vector input
        specifies a univariate mixture.
    sigmas : array_like, shape (k,) or (q, q, k)
        Component variances for vector ``means``, otherwise covariance matrices.
        These are variances, not standard deviations. Every covariance must be
        finite, symmetric and positive definite, even for zero-weight components.
    weights : array_like, shape (k,)
        Finite, non-negative weights, normalized internally. At least one must
        be positive.
    rng : int or numpy.random.Generator, optional
        Seed or generator; existing generators advance in place. ``None``
        creates a new independently seeded generator.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)`` when ``means`` is a vector, otherwise ``(n, q)``.
        Matrix-shaped univariate means therefore retain a single output column.

    Raises
    ------
    ValueError
        If dimensions, means, covariances, weights or sample size are invalid.

    See Also
    --------
    dmixnorm : Evaluate the corresponding mixture density.

    Examples
    --------
    >>> rmixnorm(5, means=[-1, 1], sigmas=[1, 2], weights=[2, 3], rng=1).shape
    (5,)
    """
    n = _validate_order("n", n)

    means, factors, univariate = _coerce_mixture_inputs(means, sigmas)
    q, k = means.shape
    weights = _normalise_weights(np.asarray(weights, dtype=float), k)
    generator = get_rng(rng)

    components = generator.choice(k, size=n, replace=True, p=weights)
    out = np.empty((n, q), dtype=float)
    for component in range(k):
        rows = components == component
        count = int(rows.sum())
        if count:
            out[rows, :] = (generator.standard_normal((count, q))
                            @ factors[component].T + means[:, component])

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


def _log_multivariate_normal_density(x: np.ndarray, mean: np.ndarray, factor: np.ndarray) -> np.ndarray:
    logdet = 2 * np.sum(np.log(np.diag(factor)))
    centered = x - mean
    solved = np.linalg.solve(factor, centered.T)
    quadratic = np.sum(solved * solved, axis=0)
    q = x.shape[1]
    return -0.5 * (q * np.log(2.0 * np.pi) + logdet + quadratic)


def _logsumexp(values: np.ndarray, axis: int) -> np.ndarray:
    return np.logaddexp.reduce(values, axis=axis)


def dmixnorm(
    x: np.ndarray | Sequence[float],
    means: np.ndarray | Sequence[float],
    sigmas: np.ndarray | Sequence[float],
    weights: Sequence[float],
    *,
    log: bool = False,
) -> np.ndarray:
    """Evaluate a Gaussian mixture density or log density.

    Parameters
    ----------
    x : array_like
        Univariate observations, one vector of length ``q``, or an ``(n, q)``
        matrix with one multivariate observation per row.
    means : array_like, shape (k,) or (q, k)
        Finite component means, using the layout accepted by ``rmixnorm``.
    sigmas : array_like, shape (k,) or (q, q, k)
        Positive variances or symmetric positive-definite covariance matrices,
        not standard deviations. Zero-weight components are also validated.
    weights : array_like, shape (k,)
        Finite, non-negative weights with at least one positive value.
        Weights are normalized internally.
    log : bool, default False
        Return log densities. Component log densities are combined stably,
        avoiding premature underflow in distribution tails.

    Returns
    -------
    numpy.ndarray
        One value per observation, normally shape ``(n,)``. A scalar ``x``
        with vector-shaped ``means`` returns a zero-dimensional array.

    Raises
    ------
    ValueError
        If observation dimensions or mixture parameters are invalid.

    See Also
    --------
    rmixnorm : Sample from the corresponding mixture.

    Examples
    --------
    >>> dmixnorm([0, 1], means=[0], sigmas=[1], weights=[1]).shape
    (2,)
    """
    means, factors, univariate = _coerce_mixture_inputs(means, sigmas)
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
                + _log_multivariate_normal_density(observations, means[:, component], factors[component])
            )
    log_density = _logsumexp(component_log_density, axis=1)
    density = log_density if log else np.exp(log_density)
    if univariate and np.asarray(x).ndim == 0:
        return np.asarray(density[0])
    return density


def _prepare_ts_mixture(n, means_ar_par_list, sigmas_list, weights):
    if not len(means_ar_par_list) or len(sigmas_list) != len(means_ar_par_list):
        raise ValueError("sigmas_list must match a nonempty means_ar_par_list")
    params = [np.asarray(values, dtype=float).reshape(-1) for values in means_ar_par_list]
    if any(values.size < 2 for values in params):
        raise ValueError("each AR parameter vector must include an intercept and at least one lag")
    if any(np.any(~np.isfinite(values)) for values in params):
        raise ValueError("AR parameters must be finite")
    sigmas = np.column_stack([np.asarray(values, dtype=float).reshape(-1) for values in sigmas_list])
    if sigmas.shape != (n, len(params)):
        raise ValueError("each sigma vector must have length n")
    if np.any(~np.isfinite(sigmas)) or np.any(sigmas <= 0):
        raise ValueError("sigmas must be positive and finite")
    lags = [np.arange(1, values.size) for values in params]
    return params, lags, np.sqrt(sigmas), _normalise_weights(np.asarray(weights), len(params))


def rmixnorm_ts(
    n: int,
    means_ar_par_list: Sequence[Sequence[float]],
    sigmas_list: Sequence[Sequence[float]],
    weights: Sequence[float],
    *,
    yinit: float = 0.0,
    rng: SeedLike = None,
) -> np.ndarray:
    """Simulate a process with Gaussian-mixture conditional distributions.

    Parameters
    ----------
    n : int
        Non-negative total path length, including the initial observations.
    means_ar_par_list : sequence of array_like
        Component vectors ``[intercept, lag_1, ..., lag_p]``. Each component
        must have at least one lag; orders may differ between components.
    sigmas_list : sequence of array_like
        One positive, finite variance vector of length ``n`` per component.
        Values are variances, not standard deviations.
    weights : array_like
        Finite, non-negative component weights, normalized internally, with
        at least one positive value.
    yinit : float, default 0
        Value assigned to the first ``min(n, max_order)`` observations. These
        are returned as part of the path; no burn-in is discarded.
    rng : int or numpy.random.Generator, optional
        Seed or generator; existing generators advance in place. ``None``
        creates a new independently seeded generator.

    Returns
    -------
    numpy.ndarray
        Floating-point path of shape ``(n,)``. At each simulated time a
        component is selected, then its mean is computed from the common past.

    Raises
    ------
    ValueError
        If parameter dimensions, variances, weights, length or ``yinit`` are
        invalid. Parameters are validated even when the path is too short to
        enter the recursion. Stationarity is not enforced.

    See Also
    --------
    dmixnorm_ts : Evaluate conditional densities along an observed path.
    mar_model : Construct a MAR model with a separate burn-in control.

    Examples
    --------
    >>> rmixnorm_ts(5, [[0, 0.5]], [[1] * 5], [1], rng=1).shape
    (5,)
    """
    n = _validate_order("n", n)
    ar_params, lags, scales, weights = _prepare_ts_mixture(
        n, means_ar_par_list, sigmas_list, weights)
    if not np.isfinite(yinit):
        raise ValueError("yinit must be finite")
    k = len(ar_params)

    generator = get_rng(rng)
    y = np.full(n, yinit, dtype=float)
    max_lag = max(len(indices) for indices in lags)
    for index in range(max_lag, n):
        component = generator.choice(k, p=weights)
        mean = np.dot(ar_params[component], np.r_[1.0, y[index - lags[component]]])
        y[index] = generator.normal(mean, scales[index, component])

    return y


def dmixnorm_ts(
    y: Sequence[float] | np.ndarray,
    means_ar_par_list: Sequence[Sequence[float]],
    sigmas_list: Sequence[Sequence[float]],
    weights: Sequence[float],
    *,
    log: bool = False,
) -> np.ndarray:
    """Evaluate conditional Gaussian-mixture densities along an observed path.

    Parameters
    ----------
    y : array_like
        Univariate observations, flattened to a one-dimensional array.
    means_ar_par_list : sequence of array_like
        Component vectors ``[intercept, lag_1, ..., lag_p]`` with at least one
        lag per component, as in ``rmixnorm_ts``.
    sigmas_list : sequence of array_like
        Positive, finite variance vectors, one per component, each of length
        ``len(y)``. These are variances, not standard deviations.
    weights : array_like
        Finite, non-negative weights, normalized internally. At least one must
        be positive; zero-weight component parameters are still validated.
    log : bool, default False
        Return conditional log densities instead of densities.

    Returns
    -------
    numpy.ndarray
        Shape ``(len(y),)``. The first ``min(len(y), max_order)`` entries are
        placeholders: one for densities, zero for log densities. Only later
        entries are evaluated conditional on the observed past. Exclude these
        placeholders when constructing a conditional likelihood.

    Raises
    ------
    ValueError
        If parameter dimensions, variances or weights are invalid.

    See Also
    --------
    rmixnorm_ts : Simulate using the same conditional-distribution convention.
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.size
    ar_params, lags, scales, weights = _prepare_ts_mixture(
        n, means_ar_par_list, sigmas_list, weights)
    max_lag = max(len(indices) for indices in lags)
    out_log = np.zeros(n, dtype=float)
    if n > max_lag:
        past = np.r_[0.0, y[:-1]]
        component_logs = []
        # Observed lags are fixed, unlike simulation's recursive state.
        for component in np.flatnonzero(weights):
            params = ar_params[component]
            means = params[0] + np.convolve(past, params[1:])[:n]
            scale = scales[max_lag:, component]
            z = (y[max_lag:] - means[max_lag:]) / scale
            component_logs.append(
                np.log(weights[component]) - np.log(scale)
                - 0.5 * (np.log(2 * np.pi) + z * z))
        out_log[max_lag:] = np.logaddexp.reduce(component_logs, axis=0)

    return out_log if log else np.exp(out_log)
