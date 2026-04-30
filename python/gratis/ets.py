"""ETS model structure and simulation for the Python gratis implementation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any
import warnings

import numpy as np

from ._random import SeedLike, get_rng


def _load_statsmodels_ets() -> tuple[type[Any] | None, Exception | None]:
    try:
        from statsmodels.tsa.exponential_smoothing.ets import ETSModel as StatsmodelsETSModel
    except Exception as exc:
        return None, exc
    return StatsmodelsETSModel, None


def _component(value: str, allowed: set[str], name: str) -> str:
    value = value.upper()
    if value not in allowed:
        raise ValueError(f"{name} must be one of {sorted(allowed)}")
    return value


def ets_admissible(
    alpha: float,
    beta: float | None = None,
    gamma: float | None = None,
    phi: float | None = None,
    m: int = 1,
) -> bool:
    """Check ETS smoothing-parameter admissibility.

    This is a direct translation of the helper used by the R implementation.
    """
    if phi is None:
        phi = 1.0
    if phi < 0 or phi > 1 + 1e-8:
        return False
    if gamma is None:
        if alpha < 1 - 1 / phi or alpha > 1 + 1 / phi:
            return False
        if beta is not None and (beta < alpha * (phi - 1) or beta > (1 + phi) * (2 - alpha)):
            return False
    elif m > 1:
        if beta is None:
            beta = 0.0
        if gamma < max(1 - 1 / phi - alpha, 0) or gamma > 1 + 1 / phi - alpha:
            return False
        if alpha < 1 - 1 / phi - gamma * (1 - m + phi + phi * m) / (2 * phi * m):
            return False
        if beta < -(1 - phi) * (gamma / m + alpha):
            return False

        polynomial = np.array(
            [
                phi * (1 - alpha - gamma),
                alpha + beta - alpha * phi + gamma - 1,
                *([alpha + beta - alpha * phi] * max(m - 2, 0)),
                alpha + beta - phi,
                1.0,
            ],
            dtype=float,
        )
        roots = np.roots(polynomial[::-1])
        if roots.size and np.max(np.abs(roots)) > 1 + 1e-10:
            return False
    return True


@dataclass(frozen=True)
class ETSModel:
    """ETS state-space model.

    The implementation supports the component set produced by the R
    `ets_model()` helper: additive or multiplicative error, no or additive
    trend, and no, additive, or multiplicative seasonality.
    """

    frequency: int = 1
    error: str = "A"
    trend: str = "N"
    seasonal: str = "N"
    damped: bool = False
    alpha: float = 0.2
    beta: float | None = None
    gamma: float | None = None
    phi: float | None = None
    level: float = 0.0
    slope: float = 0.0
    season: np.ndarray | None = None
    sigma: float = 1.0
    backend: str = "auto"

    def __post_init__(self) -> None:
        object.__setattr__(self, "error", _component(self.error, {"A", "M"}, "error"))
        object.__setattr__(self, "trend", _component(self.trend, {"N", "A"}, "trend"))
        object.__setattr__(self, "seasonal", _component(self.seasonal, {"N", "A", "M"}, "seasonal"))
        season = np.array([], dtype=float) if self.season is None else np.asarray(self.season, dtype=float).reshape(-1)
        object.__setattr__(self, "season", season)
        backend = self.backend.lower()
        if backend not in {"auto", "statsmodels", "numpy"}:
            raise ValueError("backend must be one of 'auto', 'statsmodels', or 'numpy'")
        object.__setattr__(self, "backend", backend)

    @property
    def method(self) -> str:
        trend = f"{self.trend}d" if self.trend == "A" and self.damped else self.trend
        return f"ETS({self.error},{trend},{self.seasonal})"

    def simulate(self, n: int = 100, *, rng: SeedLike = None) -> np.ndarray:
        """Simulate one sample path from the ETS model."""
        n = int(n)
        if n < 0:
            raise ValueError("n must be non-negative")
        if n == 0:
            return np.array([], dtype=float)
        if self.sigma <= 0:
            raise ValueError("sigma must be positive")
        if self.seasonal != "N" and self.frequency <= 1:
            raise ValueError("seasonal models must have frequency greater than 1")

        if self.backend in {"auto", "statsmodels"}:
            model_class, error = _load_statsmodels_ets()
            if model_class is not None:
                return self._simulate_statsmodels(model_class, n=n, rng=rng)
            if self.backend == "statsmodels":
                raise ImportError(
                    "statsmodels could not be imported for ETS simulation. "
                    "Install the Python ETS extra with `python -m pip install -e './python[ets]'`."
                ) from error

        return self._simulate_numpy(n=n, rng=rng)

    def _simulate_statsmodels(
        self,
        model_class: type[Any],
        *,
        n: int,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Simulate with statsmodels' ETSModel backend."""
        generator = get_rng(rng)
        error = "add" if self.error == "A" else "mul"
        trend = None if self.trend == "N" else "add"
        seasonal = None if self.seasonal == "N" else ("add" if self.seasonal == "A" else "mul")
        endog_value = max(float(self.level), 1.0) if "mul" in {error, seasonal} else float(self.level)
        endog = np.repeat(endog_value, max(2, self.frequency)).astype(float)

        kwargs: dict[str, Any] = {
            "error": error,
            "trend": trend,
            "damped_trend": bool(self.damped),
            "seasonal": seasonal,
            "initialization_method": "known",
            "initial_level": float(self.level),
        }
        if self.trend != "N":
            kwargs["initial_trend"] = float(self.slope)
        if self.seasonal != "N":
            kwargs["seasonal_periods"] = int(self.frequency)
            kwargs["initial_seasonal"] = np.asarray(self.season, dtype=float)

        model = model_class(endog, **kwargs)
        params = [float(self.alpha)]
        if self.trend != "N":
            params.append(0.0 if self.beta is None else float(self.beta))
        if self.seasonal != "N":
            params.append(0.0 if self.gamma is None else float(self.gamma))
        if self.damped:
            params.append(1.0 if self.phi is None else float(self.phi))

        # The statsmodels object requires an endogenous sample even when we
        # simulate from known states. That synthetic sample can produce
        # meaningless likelihood warnings while the result wrapper is built.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, module="statsmodels")
            results = model.smooth(np.asarray(params, dtype=float))
        errors = generator.normal(loc=0.0, scale=float(self.sigma), size=(n, 1))
        simulated = results.simulate(
            nsimulations=n,
            anchor="start",
            random_errors=errors,
        )
        return np.asarray(simulated, dtype=float).reshape(-1)

    def _simulate_numpy(self, n: int = 100, *, rng: SeedLike = None) -> np.ndarray:
        """Fallback ETS simulation used when statsmodels is unavailable."""
        generator = get_rng(rng)
        level = float(self.level)
        slope = float(self.slope if self.trend != "N" else 0.0)
        damping = float(self.phi if self.damped and self.phi is not None else 1.0)
        beta = 0.0 if self.beta is None else float(self.beta)
        gamma = 0.0 if self.gamma is None else float(self.gamma)
        season = np.asarray(self.season, dtype=float).copy()
        if self.seasonal != "N" and season.size != self.frequency:
            raise ValueError("season must have length equal to frequency")

        y = np.empty(n, dtype=float)
        for t in range(n):
            seasonal_state = 0.0 if self.seasonal == "N" else season[t % self.frequency]
            trend_state = damping * slope if self.trend == "A" else 0.0
            base = level + trend_state
            fitted = base
            if self.seasonal == "A":
                fitted = base + seasonal_state
            elif self.seasonal == "M":
                fitted = base * seasonal_state

            if self.error == "A":
                absolute_error = float(generator.normal(0.0, self.sigma))
                y[t] = fitted + absolute_error
                relative_error = absolute_error / fitted if abs(fitted) > 1e-12 else 0.0
            else:
                relative_error = float(generator.normal(0.0, self.sigma))
                y[t] = fitted * (1.0 + relative_error)
                absolute_error = y[t] - fitted

            if self.error == "M":
                new_level = base * (1.0 + self.alpha * relative_error)
                trend_error = base * relative_error
            elif self.seasonal == "M" and abs(seasonal_state) > 1e-12:
                new_level = base + self.alpha * absolute_error / seasonal_state
                trend_error = absolute_error / seasonal_state
            else:
                new_level = base + self.alpha * absolute_error
                trend_error = absolute_error

            if self.trend == "A":
                new_slope = trend_state + beta * trend_error
            else:
                new_slope = 0.0

            if self.seasonal == "A":
                season[t % self.frequency] = seasonal_state + gamma * absolute_error
            elif self.seasonal == "M":
                if self.error == "M":
                    season[t % self.frequency] = seasonal_state * (1.0 + gamma * relative_error)
                elif abs(base) > 1e-12:
                    season[t % self.frequency] = seasonal_state + gamma * absolute_error / base

            level = new_level
            slope = new_slope

        return y


def ets_model(
    frequency: int = 1,
    error: str | None = None,
    trend: str | None = None,
    seasonal: str | None = None,
    alpha: float | None = None,
    beta: float | None = None,
    gamma: float | None = None,
    phi: float | None = None,
    level: float | None = None,
    slope: float | None = None,
    season: Sequence[float] | np.ndarray | None = None,
    damped: bool | None = None,
    sigma: float | None = None,
    backend: str = "auto",
    *,
    rng: SeedLike = None,
) -> ETSModel:
    """Specify an ETS state-space model and return a simulatable object."""
    generator = get_rng(rng)
    frequency = int(frequency)
    if frequency < 1:
        raise ValueError("frequency must be a positive integer")

    if error is None:
        error = str(generator.choice(["A", "M"]))
    error = _component(error, {"A", "M"}, "error")

    if sigma is None:
        sigma2 = float(generator.uniform(1.0, 5.0)) if error == "A" else float(generator.uniform(1e-4, 0.05))
        sigma = float(np.sqrt(sigma2))
    elif sigma <= 0:
        raise ValueError("sigma must be positive")

    if trend is None:
        trend = str(generator.choice(["N", "A"]))
    trend = _component(trend, {"N", "A"}, "trend")

    if damped is None:
        damped = bool(generator.choice([True, False])) if trend == "A" else False
    if trend == "N":
        damped = False

    if seasonal is None:
        seasonal = str(generator.choice(["N", "A", "M"])) if frequency > 1 else "N"
    seasonal = _component(seasonal, {"N", "A", "M"}, "seasonal")
    if seasonal != "N" and frequency <= 1:
        raise ValueError("seasonal models must have frequency greater than 1")
    m = frequency if seasonal != "N" else 1

    supplied = {
        "alpha": alpha is not None,
        "beta": beta is not None,
        "gamma": gamma is not None,
        "phi": phi is not None,
    }
    for _ in range(10_000):
        alpha_candidate = float(alpha) if supplied["alpha"] else float(generator.uniform(0.0, 1.0))
        beta_candidate = None
        if trend != "N":
            beta_candidate = float(beta) if supplied["beta"] else float(generator.uniform(0.0, alpha_candidate))
        gamma_candidate = None
        if seasonal != "N":
            gamma_candidate = float(gamma) if supplied["gamma"] else float(generator.uniform(0.0, 1.0 - alpha_candidate))
        phi_candidate = None
        if trend == "A" and damped:
            phi_candidate = float(phi) if supplied["phi"] else float(generator.uniform(0.8, 0.98))
        if ets_admissible(alpha_candidate, beta_candidate, gamma_candidate, phi_candidate, m):
            alpha, beta, gamma, phi = alpha_candidate, beta_candidate, gamma_candidate, phi_candidate
            break
    else:
        raise ValueError("could not find admissible ETS smoothing parameters")

    if level is None:
        level = float(generator.uniform(-1.0, 1.0)) if error == "A" else float(generator.uniform(2.0, 10.0))
    if trend != "N" and slope is None:
        slope = float(generator.uniform(-1.0, 1.0))
    elif trend == "N":
        slope = 0.0

    if season is None and seasonal == "A":
        first = generator.uniform(-1.0, 1.0, size=m - 1)
        season = np.concatenate((first, np.array([-float(first.sum())])))
    elif season is None and seasonal == "M":
        first = generator.uniform(0.5, 1.5, size=m - 1)
        season = np.concatenate((first, np.array([m - float(first.sum())])))
    elif season is None:
        season = np.array([], dtype=float)
    else:
        season = np.asarray(season, dtype=float).reshape(-1)
        if seasonal != "N" and season.size != m:
            raise ValueError("season must have length equal to frequency")

    return ETSModel(
        frequency=frequency,
        error=error,
        trend=trend,
        seasonal=seasonal,
        damped=bool(damped),
        alpha=float(alpha),
        beta=None if beta is None else float(beta),
        gamma=None if gamma is None else float(gamma),
        phi=None if phi is None else float(phi),
        level=float(level),
        slope=float(slope),
        season=np.asarray(season, dtype=float),
        sigma=float(sigma),
        backend=backend,
    )
