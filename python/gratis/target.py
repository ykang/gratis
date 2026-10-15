"""Target-feature generation workflows."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import warnings

from ._random import SeedLike, get_rng
from .arima import _validate_order
from .mar import mar_model


def scalets(x: Sequence[float] | np.ndarray) -> np.ndarray:
    """Center and sample-standardize finite data without squaring large values.

    Parameters
    ----------
    x : array_like
        Finite observations, flattened to a one-dimensional array. Empty and
        exactly constant inputs, including singletons, are returned unchanged
        in value. The input is not modified.

    Returns
    -------
    numpy.ndarray
        Floating-point vector. Nonconstant data have mean zero and sample
        standard deviation one, up to numerical precision; ``ddof=1`` is used.

    Raises
    ------
    ValueError
        If any observation is nonfinite.

    Examples
    --------
    >>> scalets([1000001, 1000002, 1000003]).tolist()
    [-1.0, 0.0, 1.0]
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size == 0:
        return x.copy()
    if np.any(~np.isfinite(x)):
        raise ValueError("x must contain only finite values")
    if np.all(x == x[0]):
        return x.copy()
    # Shift before rescaling to retain small variations around large offsets.
    with np.errstate(over="ignore"):
        shifted = x - x[0]
    if np.any(~np.isfinite(shifted)):
        bounded = x / np.max(np.abs(x))
        shifted = bounded - bounded[0]
    shifted /= np.max(np.abs(shifted))
    centered = shifted - np.mean(shifted)
    return centered / np.sqrt(np.sum(centered * centered) / (x.size - 1))


def _reshape_column_major(values: np.ndarray, rows: int, columns: int) -> np.ndarray:
    if rows == 0:
        return np.zeros((0, columns), dtype=float)
    return np.asarray(values, dtype=float).reshape((rows, columns), order="F")


def _params_to_mar(
    params: np.ndarray,
    seasonal_periods: np.ndarray,
    k: int,
) -> object:
    p = 3
    P = 2 if np.max(seasonal_periods) > 1 else 0

    d = np.rint(params[:k]).astype(int)
    D = np.rint(params[k : 2 * k]).astype(int)
    phi_start = 2 * k
    phi = _reshape_column_major(params[phi_start : phi_start + p * k], p, k)
    Phi_start = phi_start + p * k
    Phi = _reshape_column_major(params[Phi_start : Phi_start + P * k], P, k)
    constants_start = Phi_start + P * k
    constants = params[constants_start : constants_start + k]
    sigmas_start = constants_start + k
    sigmas = np.maximum(params[sigmas_start : sigmas_start + k], 1e-8)
    weights_start = sigmas_start + k
    raw_weights = params[weights_start : weights_start + k - 1]
    weights = np.diff(np.sort(np.r_[0.0, raw_weights, 1.0]))

    return mar_model(
        k=k,
        d=np.clip(d, 0, 2),
        D=np.clip(D, 0, 1),
        p=p,
        P=P,
        phi=phi,
        Phi=Phi,
        constants=constants,
        sigmas=sigmas,
        weights=weights,
        seasonal_periods=seasonal_periods,
    )


def _fitness_mar(
    params: np.ndarray,
    *,
    length: int,
    seasonal_periods: np.ndarray,
    k: int,
    feature_function: Callable[[np.ndarray], Sequence[float] | np.ndarray],
    target: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    model = _params_to_mar(params, seasonal_periods, k)
    with np.errstate(over="ignore", invalid="ignore"):
        series = model.simulate(n=length, rng=rng)
    if not np.all(np.isfinite(series)):
        return -np.inf, np.full(length, np.nan)
    series = scalets(series)
    features = np.asarray(feature_function(series), dtype=float).reshape(-1)
    if features.shape != target.shape:
        raise ValueError("feature_function must return one value per target feature")
    if np.any(~np.isfinite(features)):
        return -np.inf, series
    return -float(np.mean(np.abs(features - target))), series


def _target_bounds(seasonal_periods: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    p = 3
    P = 2 if np.max(seasonal_periods) > 1 else 0
    par_length = (2 + p + P + 3) * k - 1
    lower = np.zeros(par_length, dtype=float)
    upper = np.ones(par_length, dtype=float)
    upper[:k] = 2.0
    coefficient_slice = slice(2 * k, 2 * k + p * k + P * k)
    lower[coefficient_slice] = -1.5
    upper[coefficient_slice] = 1.5
    constants_slice = slice((2 + p + P) * k, (3 + p + P) * k)
    lower[constants_slice] = -3.0
    upper[constants_slice] = 3.0
    sigmas_slice = slice((3 + p + P) * k, (4 + p + P) * k)
    upper[sigmas_slice] = 5.0
    return lower, upper


def simulate_target(
    length: int = 100,
    seasonal_periods: int | Sequence[int] = 1,
    feature_function: Callable[[np.ndarray], Sequence[float] | np.ndarray] | None = None,
    target: Sequence[float] | np.ndarray | None = None,
    k: int | None = None,
    tolerance: float = 0.05,
    trace: bool = False,
    parallel: bool = False,
    *,
    rng: SeedLike = None,
    population_size: int = 100,
    max_iter: int = 200,
) -> np.ndarray:
    """Generate one series whose features are close to a target vector.

    The R package uses a custom GA object. The Python port implements the same
    search idea with a compact evolutionary loop over MAR parameters.
    If tolerance is not reached, the best series is returned with a warning
    containing its mean absolute feature error. Callback exceptions propagate;
    nonfinite features or simulated paths are rejected as invalid candidates.
    With parallel=True, candidates run in threads. Feature callbacks must be
    thread-safe; speedups depend on their releasing the Python GIL. Independent
    candidate random streams give identical serial and parallel seeded results.

    Parameters
    ----------
    length : int, default 100
        Positive number of observations in each candidate and the returned path.
    seasonal_periods : int or array_like, default 1
        Positive integer seasonal period, or one period per MAR component.
        A scalar is applied to every component.
    feature_function : callable
        Function accepting a one-dimensional candidate series and returning one
        numeric value per target feature. A scalar result is allowed for one
        feature. Candidates are sample-standardized first; targets should not
        rely on their original scale. Exactly constant candidates are unchanged.
    target : array_like
        Nonempty finite feature vector, in the callback's output order.
    k : int, optional
        Positive component count. Defaults to three for a scalar seasonal period
        or to the number of supplied periods otherwise. A period vector must
        have exactly ``k`` entries.
    tolerance : float, default 0.05
        Non-negative bound on the mean absolute difference between features
        and target. Reaching it stops the search; success is not guaranteed.
    trace : bool, default False
        Print the iteration number and best score. Scores are negative mean
        absolute feature errors, so larger values are better.
    parallel : bool, default False
        Evaluate candidates in threads. Callbacks must be thread-safe. Speedups
        depend on callback cost and whether it releases the Python GIL.
    rng : int or numpy.random.Generator, optional
        Seed or generator for the search and candidate streams. An existing
        generator advances in place; ``None`` creates a new one. Serial and
        threaded results match for a fixed seed and deterministic callback,
        but Python and R searches need not return the same path.
    population_size : int, default 100
        Number of candidates per generation; at least two.
    max_iter : int, default 200
        Positive maximum number of generations.

    Returns
    -------
    numpy.ndarray
        Best finite candidate found, shape ``(length,)``, on the scale used for
        feature evaluation. No fitted MAR model or convergence object is returned.

    Raises
    ------
    ValueError
        If arguments are invalid or the callback returns an incompatible length.
    RuntimeError
        If the search cannot produce a finite, valid candidate.

    Warns
    -----
    RuntimeWarning
        If the returned path misses ``tolerance``. The warning reports its mean
        absolute feature error. Exceptions raised by callbacks are not suppressed.

    See Also
    --------
    generate_target : Repeat the search to obtain several paths.
    scalets : Standardization applied before feature evaluation.

    Examples
    --------
    >>> y = simulate_target(length=20, k=1, feature_function=np.mean,
    ...                     target=[0], tolerance=1e-8, population_size=4,
    ...                     max_iter=1, rng=1)
    >>> y.shape
    (20,)
    """
    if feature_function is None:
        raise ValueError("feature_function must be provided")
    if target is None:
        raise ValueError("target must be provided")
    if not callable(feature_function):
        raise ValueError("feature_function must be callable")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be non-negative and finite")
    length = _validate_order("length", length)
    if length < 1:
        raise ValueError("length must be positive")
    population_size = _validate_order("population_size", population_size)
    max_iter = _validate_order("max_iter", max_iter)
    if population_size < 2:
        raise ValueError("population_size must be at least 2")
    if max_iter < 1:
        raise ValueError("max_iter must be positive")

    generator = get_rng(rng)
    seasonal_periods_arr = np.asarray(seasonal_periods, dtype=float).reshape(-1)
    if (seasonal_periods_arr.size == 0 or np.any(~np.isfinite(seasonal_periods_arr))
            or np.any(seasonal_periods_arr < 1)
            or np.any(seasonal_periods_arr != np.floor(seasonal_periods_arr))):
        raise ValueError("seasonal_periods must contain positive integers")
    seasonal_periods_arr = seasonal_periods_arr.astype(int)
    if k is None:
        k = 3 if seasonal_periods_arr.size == 1 else seasonal_periods_arr.size
    k = _validate_order("k", k)
    if k < 1:
        raise ValueError("k must be positive")
    target_arr = np.asarray(target, dtype=float).reshape(-1)
    if target_arr.size == 0 or np.any(~np.isfinite(target_arr)):
        raise ValueError("target must be a nonempty finite numeric vector")

    lower, upper = _target_bounds(seasonal_periods_arr, k)
    span = upper - lower
    population = generator.uniform(lower, upper, size=(population_size, lower.size))
    best_score = -np.inf
    best_series = np.full(length, np.nan)

    def evaluate(item):
        candidate, seed = item
        return _fitness_mar(
            candidate,
            length=length,
            seasonal_periods=seasonal_periods_arr,
            k=k,
            feature_function=feature_function,
            target=target_arr,
            rng=get_rng(int(seed)),
        )

    elite_count = max(1, population_size // 5)
    # No threads are started by the executor unless parallel mapping is used.
    with ThreadPoolExecutor() as executor:
        for iteration in range(max_iter):
            seeds = generator.integers(0, np.iinfo(np.int64).max, size=population_size)
            items = zip(population, seeds)
            results = list(executor.map(evaluate, items) if parallel else map(evaluate, items))
            scores = np.asarray([result[0] for result in results])
            if not np.any(np.isfinite(scores)):
                if np.isfinite(best_score):
                    break
                raise RuntimeError("target search did not produce a finite series")
            current = int(np.nanargmax(scores))
            if scores[current] > best_score:
                best_score, best_series = results[current]
            if trace:
                print(f"iteration={iteration + 1} score={best_score:.6g}")
            if -best_score <= tolerance:
                break

            elite_indices = np.argsort(scores)[-elite_count:]
            elites = population[elite_indices]
            next_population = np.empty_like(population)
            next_population[:elite_count] = elites
            for row in range(elite_count, population_size):
                parent_a, parent_b = elites[generator.integers(0, elite_count, size=2)]
                child = 0.5 * (parent_a + parent_b)
                mutation = generator.normal(0.0, 0.1, size=lower.size) * span
                reset = generator.random(lower.size) < 0.05
                child = child + mutation
                child[reset] = generator.uniform(lower[reset], upper[reset])
                next_population[row] = np.clip(child, lower, upper)
            population = next_population

    if not np.isfinite(best_score) or not np.all(np.isfinite(best_series)):
        raise RuntimeError("target search did not produce a finite series")
    if -best_score > tolerance:
        warnings.warn(
            f"target tolerance not reached: mean absolute feature error={-best_score:.6g}, "
            f"tolerance={tolerance:.6g}",
            RuntimeWarning,
            stacklevel=2,
        )
    return np.asarray(best_series, dtype=float)


def generate_target(
    length: int = 100,
    nseries: int = 10,
    seasonal_periods: int | Sequence[int] = 1,
    feature_function: Callable[[np.ndarray], Sequence[float] | np.ndarray] | None = None,
    target: Sequence[float] | np.ndarray | None = None,
    k: int | None = None,
    tolerance: float = 0.05,
    trace: bool = False,
    parallel: bool = False,
    *,
    rng: SeedLike = None,
    population_size: int = 100,
    max_iter: int = 200,
) -> np.ndarray:
    """Repeat target-feature searches, returning one path per column.

    Parameters
    ----------
    length : int, default 100
        Positive number of observations in each returned path.
    nseries : int, default 10
        Positive number of separate searches to perform.
    seasonal_periods : int or array_like, default 1
        Seasonal period or per-component periods, as in ``simulate_target``.
    feature_function : callable
        Feature callback applied to sample-standardized candidates. Must return
        one numeric value per target feature.
    target : array_like
        Nonempty finite target vector in callback order.
    k : int, optional
        MAR component count, with the defaults used by ``simulate_target``.
    tolerance : float, default 0.05
        Non-negative mean absolute feature-error tolerance for each search.
    trace, parallel : bool, default False
        Enable progress output or threaded candidate evaluation, respectively.
        Searches themselves are performed sequentially.
    rng : int or numpy.random.Generator, optional
        Seed or shared generator. Initialized once per batch; successive
        searches use successive draws rather than restarting the same seed.
        ``None`` creates a new independently seeded generator.
    population_size : int, default 100
        Number of candidates per generation, at least two.
    max_iter : int, default 200
        Positive generation limit per search.

    Returns
    -------
    numpy.ndarray
        Shape ``(length, nseries)``, including for a single series. Each column
        is the best path from its own search, not necessarily within tolerance.

    Notes
    -----
    Tolerance-miss warnings and callback errors follow ``simulate_target``.
    Any failed search propagates its error; a partial batch is not returned.

    See Also
    --------
    simulate_target : Search objective, validation and convergence semantics.
    """
    nseries = _validate_order("nseries", nseries)
    if nseries < 1:
        raise ValueError("nseries must be positive")
    generator = get_rng(rng)
    out = np.empty((int(length), nseries), dtype=float)
    for index in range(nseries):
        out[:, index] = simulate_target(
            length=length,
            seasonal_periods=seasonal_periods,
            feature_function=feature_function,
            target=target,
            k=k,
            tolerance=tolerance,
            trace=trace,
            parallel=parallel,
            rng=generator,
            population_size=population_size,
            max_iter=max_iter,
        )
    return out
