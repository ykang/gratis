"""ETS model structure and simulation for the Python gratis implementation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache, partial
from typing import Any
import warnings

import numpy as np

from ._random import SeedLike, get_rng
from .arima import _validate_order, _validate_sigma


def _component(value: str, allowed: set[str], name: str) -> str:
    value = value.upper()
    if value not in allowed:
        raise ValueError(f"{name} must be one of {sorted(allowed)}")
    return value


@lru_cache(maxsize=128)
def ets_admissible(
    alpha: float,
    beta: float | None = None,
    gamma: float | None = None,
    phi: float | None = None,
    m: int = 1,
) -> bool:
    """Check ETS smoothing-parameter admissibility.

    Parameters
    ----------
    alpha : float
        Level smoothing parameter.
    beta, gamma : float, optional
        Trend and seasonal smoothing parameters. ``None`` omits that component.
    phi : float, optional
        Damping parameter; ``None`` is treated as one.
    m : int, default 1
        Seasonal period used when a seasonal smoothing parameter is supplied.

    Returns
    -------
    bool
        Whether the parameters pass the R-compatible ETS admissibility test.
        This low-level check does not replace the usual smoothing bounds or
        model-dimension validation performed by ``ets_model`` and ``ETSModel``.
    """
    if phi is None:
        phi = 1.0
    if any(not np.isfinite(x) for x in (alpha, beta, gamma, phi) if x is not None):
        return False
    if phi <= 0 or phi > 1 + 1e-8:
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


def _validate_smoothing(alpha, beta, gamma, phi):
    for name, value in (("alpha", alpha), ("beta", beta), ("gamma", gamma), ("phi", phi)):
        if value is None:
            continue
        if np.ndim(value) != 0 or not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError(f"{name} must be finite and between 0 and 1")
        if name in {"alpha", "phi"} and value == 0:
            raise ValueError(f"{name} must be positive")
        if name == "alpha" and value == 1:
            raise ValueError("alpha must be less than 1")


@dataclass(frozen=True)
class ETSModel:
    """Store a specified ETS innovations state-space model.

    Parameters
    ----------
    frequency : int, default 1
        Positive seasonal period. Seasonal models require a value above one.
    error : {"A", "M"}, default "A"
        Additive or multiplicative observation errors.
    trend : {"N", "A"}, default "N"
        No trend or additive trend.
    seasonal : {"N", "A", "M"}, default "N"
        No, additive or multiplicative seasonality.
    damped : bool, default False
        Whether an additive trend is damped. Disabled when there is no trend.
    alpha : float, default 0.2
        Level smoothing parameter, strictly between zero and one.
    beta, gamma, phi : float, optional
        Trend smoothing, seasonal smoothing and damping parameters. Active
        missing smoothing parameters default to zero; active missing damping
        defaults to one. Inactive parameters are ignored.
    level, slope : float, default 0
        Finite initial level and additive trend state. Multiplicative models
        require a positive level. Slope is unused when there is no trend.
    season : array_like, optional
        Finite initial seasonal states, oldest to newest, with length equal to
        ``frequency`` for seasonal models. Multiplicative factors must be positive.
    sigma : float, default 1
        Innovation standard deviation; dimensionless for multiplicative error.
        Its square must lie in the positive normal floating-point range.
    backend : {"auto", "statsmodels"}, default "auto"
        Both choices use statsmodels; there is no NumPy fallback.

    Notes
    -----
    Parameters are validated at construction. Use ``ets_model`` to randomly
    sample missing components, smoothing parameters or initial states instead.
    Seasonal states are copied, but their array remains writable.
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
        frequency = _validate_order("frequency", self.frequency)
        if frequency < 1:
            raise ValueError("frequency must be a positive integer")
        object.__setattr__(self, "frequency", frequency)
        object.__setattr__(self, "error", _component(self.error, {"A", "M"}, "error"))
        object.__setattr__(self, "trend", _component(self.trend, {"N", "A"}, "trend"))
        object.__setattr__(self, "seasonal", _component(self.seasonal, {"N", "A", "M"}, "seasonal"))
        season = np.array([], dtype=float) if self.season is None else np.asarray(self.season, dtype=float).reshape(-1)
        if self.seasonal != "N" and (frequency <= 1 or season.size != frequency):
            raise ValueError("seasonal models require frequency greater than 1 and season of that length")
        if np.any(~np.isfinite(season)):
            raise ValueError("season must be finite")
        if self.seasonal == "M" and np.any(season <= 0):
            raise ValueError("multiplicative seasonal factors must be positive")
        object.__setattr__(self, "season", season.copy())
        for name in ("level", "slope", "sigma"):
            value = getattr(self, name)
            if np.ndim(value) != 0 or not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
        _validate_sigma(self.sigma)
        if (self.error == "M" or self.seasonal == "M") and self.level <= 0:
            raise ValueError("multiplicative ETS models require a positive initial level")
        if not isinstance(self.damped, (bool, np.bool_)):
            raise ValueError("damped must be a boolean")
        if self.trend == "N":
            object.__setattr__(self, "damped", False)
        beta = None if self.trend == "N" else (0.0 if self.beta is None else self.beta)
        gamma = None if self.seasonal == "N" else (0.0 if self.gamma is None else self.gamma)
        phi = None if not self.damped else (1.0 if self.phi is None else self.phi)
        _validate_smoothing(self.alpha, beta, gamma, phi)
        if self.alpha is None:
            raise ValueError("alpha must be provided")
        if (beta is not None and beta > self.alpha) or (gamma is not None and gamma > 1 - self.alpha):
            raise ValueError("beta must be <= alpha and gamma must be <= 1 - alpha")
        alpha = float(self.alpha)
        beta = None if beta is None else float(beta)
        gamma = None if gamma is None else float(gamma)
        phi = None if phi is None else float(phi)
        if not ets_admissible(alpha, beta, gamma, phi, frequency if gamma is not None else 1):
            raise ValueError("inadmissible ETS smoothing parameters")
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "gamma", gamma)
        object.__setattr__(self, "phi", phi)
        backend = self.backend.lower()
        if backend not in {"auto", "statsmodels"}:
            raise ValueError("ETS uses statsmodels; backend must be 'auto' or 'statsmodels'")
        object.__setattr__(self, "backend", "statsmodels")

    @property
    def method(self) -> str:
        """ETS component label, such as ``ETS(A,Ad,N)`` for a damped trend."""
        trend = f"{self.trend}d" if self.trend == "A" and self.damped else self.trend
        return f"ETS({self.error},{trend},{self.seasonal})"

    def simulate(self, n: int = 100, *, rng: SeedLike = None) -> np.ndarray:
        """Simulate a path from the specified initial states using statsmodels.

        Parameters
        ----------
        n : int, default 100
            Non-negative number of output observations. No burn-in is applied.
        rng : int or numpy.random.Generator, optional
            Seed or generator for innovations. Existing generators advance in
            place; ``None`` creates a new independently seeded generator.

        Returns
        -------
        numpy.ndarray
            Floating-point path of shape ``(n,)``, or an empty array for zero.

        Raises
        ------
        ValueError
            If ``n`` or the innovation scale is invalid.
        ImportError
            If statsmodels is unavailable.

        Notes
        -----
        Multiplicative innovations are clipped below at -1 to match the R
        forecast simulator. They therefore are not unbounded Gaussian errors.
        Admissible smoothing parameters and positive initial states do not
        guarantee that every simulated path stays positive or finite.
        """
        return self._simulate(n, rng=rng)

    def _prepare_simulator(self):
        return partial(self._simulate, prepared=self._prepare_statsmodels())

    def _simulate(self, n=100, *, rng=None, prepared=None):
        n = _validate_order("n", n)
        if n == 0:
            return np.array([], dtype=float)
        _validate_sigma(self.sigma)
        if self.seasonal != "N" and self.frequency <= 1:
            raise ValueError("seasonal models must have frequency greater than 1")

        return self._simulate_statsmodels(n=n, rng=rng, prepared=prepared)

    def _prepare_statsmodels(self):
        from statsmodels.tsa.exponential_smoothing.ets import ETSModel as StatsmodelsETSModel

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

        model = StatsmodelsETSModel(endog, **kwargs)
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
        return results

    def _simulate_statsmodels(self, *, n, rng=None, prepared=None):
        """Simulate with batch-local statsmodels ETS preparation."""
        results = self._prepare_statsmodels() if prepared is None else prepared
        generator = get_rng(rng)
        errors = generator.normal(loc=0.0, scale=float(self.sigma), size=(n, 1))
        if self.error == "M":
            errors = np.maximum(errors, -1.0)
        simulated = results.simulate(
            nsimulations=n,
            anchor="start",
            random_errors=errors,
        )
        return np.asarray(simulated, dtype=float).reshape(-1)


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
    """Specify ETS using statsmodels, with R-compatible states and parameters.

    Seasonal states run from oldest to newest. Multiplicative seasonal factors
    must be positive; randomly sampled factors are normalized to mean one.
    Smoothing parameters obey the usual bounds and ETS admissibility conditions.

    Parameters
    ----------
    frequency : int, default 1
        Positive seasonal period, for example four for quarterly data.
    error : {"A", "M"}, optional
        Additive or multiplicative error, sampled when omitted.
    trend : {"N", "A"}, optional
        No trend or additive trend, sampled when omitted.
    seasonal : {"N", "A", "M"}, optional
        Seasonal component. Sampled from all three choices for periods above
        one; otherwise defaults to no seasonality.
    alpha, beta, gamma : float, optional
        Level, trend and seasonal smoothing parameters. Missing active values
        are resampled until admissible, subject to ``0 < alpha < 1``,
        ``0 <= beta <= alpha`` and ``0 <= gamma <= 1 - alpha``.
    phi : float, optional
        Damping in ``(0, 1]``, sampled uniformly from ``(0.8, 0.98)`` when active
        and omitted. Ignored when there is no damped trend.
    level, slope : float, optional
        Initial level and additive trend state. Missing values are sampled;
        multiplicative models require a positive level. No-trend models use
        zero slope.
    season : array_like, optional
        Initial seasonal states, oldest to newest, of length ``frequency``.
        Sampled additive states sum to zero; sampled multiplicative states are
        positive and have mean one. Supplied states are not renormalized.
    damped : bool, optional
        Whether to damp an additive trend. Randomly selected for additive
        trends when omitted, otherwise disabled when there is no trend.
    sigma : float, optional
        Innovation standard deviation, not variance. When omitted, the variance
        is uniform on ``(1, 5)`` for additive errors or ``(0.0001, 0.05)`` for
        multiplicative errors. Its square must lie in the positive normal
        floating-point range.
    backend : {"auto", "statsmodels"}, default "auto"
        Both select statsmodels' ETS implementation; missing dependencies raise
        an error rather than selecting another simulator.
    rng : int or numpy.random.Generator, optional
        Seed or generator for sampling missing parameters, not a stored seed
        for subsequent simulation. Existing generators advance in place;
        ``None`` creates a new independently seeded generator.

    Returns
    -------
    ETSModel
        Validated model with known initial states, ready to simulate.

    Raises
    ------
    ValueError
        If fixed parameters or initial states are invalid, or no admissible
        smoothing parameters are found within 10,000 attempts.

    See Also
    --------
    ETSModel.simulate : Simulation and multiplicative-innovation clipping.
    ets_admissible : Low-level admissibility test.

    Examples
    --------
    >>> model = ets_model(error="A", trend="N", seasonal="N", alpha=0.2,
    ...                   level=0, sigma=1)
    >>> model.simulate(12, rng=1).shape
    (12,)
    """
    generator = get_rng(rng)
    if not np.isfinite(frequency) or frequency != int(frequency):
        raise ValueError("frequency must be a positive integer")
    frequency = int(frequency)
    if frequency < 1:
        raise ValueError("frequency must be a positive integer")

    if error is None:
        error = str(generator.choice(["A", "M"]))
    error = _component(error, {"A", "M"}, "error")

    if sigma is None:
        sigma2 = float(generator.uniform(1.0, 5.0)) if error == "A" else float(generator.uniform(1e-4, 0.05))
        sigma = float(np.sqrt(sigma2))
    _validate_sigma(sigma)

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
    if trend == "N":
        beta = None
    if seasonal == "N":
        gamma = None
    if not damped:
        phi = None
    _validate_smoothing(alpha, beta, gamma, phi)
    lower = 0.0 if beta is None else beta
    upper = 1.0 if gamma is None else 1.0 - gamma
    if lower > upper or upper <= 0 or lower >= 1:
        raise ValueError("beta and gamma leave no admissible alpha")
    if alpha is not None and not lower <= alpha <= upper:
        raise ValueError("beta must be <= alpha and gamma must be <= 1 - alpha")
    all_fixed = (alpha is not None and (trend == "N" or beta is not None)
                 and (seasonal == "N" or gamma is not None)
                 and (not damped or phi is not None))
    for _ in range(10_000):
        alpha_candidate = float(alpha) if alpha is not None else float(generator.uniform(lower, upper))
        beta_candidate = None
        if trend != "N":
            beta_candidate = float(beta) if beta is not None else float(generator.uniform(0.0, alpha_candidate))
        gamma_candidate = None
        if seasonal != "N":
            gamma_candidate = float(gamma) if gamma is not None else float(generator.uniform(0.0, 1.0 - alpha_candidate))
        phi_candidate = None
        if trend == "A" and damped:
            phi_candidate = float(phi) if phi is not None else float(generator.uniform(0.8, 0.98))
        if ets_admissible(alpha_candidate, beta_candidate, gamma_candidate, phi_candidate, m):
            alpha, beta, gamma, phi = alpha_candidate, beta_candidate, gamma_candidate, phi_candidate
            break
        if all_fixed:
            raise ValueError("inadmissible ETS smoothing parameters")
    else:
        raise ValueError("could not find admissible ETS smoothing parameters")

    if level is None:
        level = float(generator.uniform(2.0, 10.0)) if error == "M" or seasonal == "M" else float(generator.uniform(-1.0, 1.0))
    if not np.isfinite(level):
        raise ValueError("level must be finite")
    if (error == "M" or seasonal == "M") and level <= 0:
        raise ValueError("multiplicative ETS models require a positive initial level")
    if trend != "N" and slope is None:
        slope = float(generator.uniform(-1.0, 1.0))
    elif trend == "N":
        slope = 0.0
    if not np.isfinite(slope):
        raise ValueError("slope must be finite")

    if season is None and seasonal == "A":
        first = generator.uniform(-1.0, 1.0, size=m - 1)
        season = np.concatenate((first, np.array([-float(first.sum())])))
    elif season is None and seasonal == "M":
        season = generator.uniform(0.5, 1.5, size=m)
        season /= season.mean()
    elif season is None:
        season = np.array([], dtype=float)
    else:
        season = np.asarray(season, dtype=float).reshape(-1)
        if seasonal != "N" and season.size != m:
            raise ValueError("season must have length equal to frequency")
    if np.any(~np.isfinite(season)):
        raise ValueError("season must be finite")
    if seasonal == "M" and np.any(season <= 0):
        raise ValueError("multiplicative seasonal factors must be positive")

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
