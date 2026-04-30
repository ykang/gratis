"""Target-feature generation workflows."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from ._random import SeedLike, get_rng
from .mar import mar_model


def scalets(x: Sequence[float] | np.ndarray) -> np.ndarray:
    """Scale a time series unless it is constant."""
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size == 0:
        return x.copy()
    std = float(np.std(x, ddof=1)) if x.size > 1 else 0.0
    if std == 0.0 or np.allclose(x, x[0]):
        return x.copy()
    return (x - float(np.mean(x))) / std


def _softmax(values: np.ndarray) -> np.ndarray:
    values = values - np.max(values)
    exp_values = np.exp(values)
    return exp_values / np.sum(exp_values)


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
    weights = _softmax(np.concatenate((raw_weights, np.array([1.0 - float(np.sum(raw_weights))]))))

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
    try:
        model = _params_to_mar(params, seasonal_periods, k)
        series = model.simulate(n=length, rng=rng)
        series = scalets(series)
        features = np.asarray(feature_function(series), dtype=float).reshape(-1)
        if features.shape != target.shape or np.any(~np.isfinite(features)):
            return -np.inf, series
        return -float(np.mean(np.abs(features - target))), series
    except Exception:
        return -np.inf, np.full(length, np.nan)


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
    """
    del parallel
    if feature_function is None:
        raise ValueError("feature_function must be provided")
    if target is None:
        raise ValueError("target must be provided")
    length = int(length)
    if length < 1:
        raise ValueError("length must be positive")
    population_size = int(population_size)
    max_iter = int(max_iter)
    if population_size < 2:
        raise ValueError("population_size must be at least 2")
    if max_iter < 1:
        raise ValueError("max_iter must be positive")

    generator = get_rng(rng)
    seasonal_periods_arr = np.asarray(seasonal_periods, dtype=int).reshape(-1)
    if seasonal_periods_arr.size == 0 or np.any(seasonal_periods_arr < 1):
        raise ValueError("seasonal_periods must contain positive integers")
    if k is None:
        k = 3 if seasonal_periods_arr.size == 1 else seasonal_periods_arr.size
    k = int(k)
    if k < 1:
        raise ValueError("k must be positive")
    target_arr = np.asarray(target, dtype=float).reshape(-1)

    lower, upper = _target_bounds(seasonal_periods_arr, k)
    span = upper - lower
    population = generator.uniform(lower, upper, size=(population_size, lower.size))
    best_score = -np.inf
    best_series = np.full(length, np.nan)

    elite_count = max(1, population_size // 5)
    for iteration in range(max_iter):
        scores = np.empty(population_size, dtype=float)
        series_list: list[np.ndarray] = []
        for row, candidate in enumerate(population):
            scores[row], series = _fitness_mar(
                candidate,
                length=length,
                seasonal_periods=seasonal_periods_arr,
                k=k,
                feature_function=feature_function,
                target=target_arr,
                rng=generator,
            )
            series_list.append(series)

        current = int(np.nanargmax(scores))
        if scores[current] > best_score:
            best_score = float(scores[current])
            best_series = series_list[current]
        if trace:
            print(f"iteration={iteration + 1} score={best_score:.6g}")
        if np.isfinite(best_score) and -best_score <= tolerance:
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

    if not np.any(np.isfinite(best_series)):
        raise RuntimeError("target search did not produce a finite series")
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
    """Generate several target-feature series as a length x nseries array."""
    nseries = int(nseries)
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
