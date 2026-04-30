"""ARIMA utilities and simulation models for the Python gratis implementation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._random import SeedLike, get_rng


def _load_statsmodels_sarimax() -> tuple[type[Any] | None, Exception | None]:
    try:
        from statsmodels.tsa.statespace.sarimax import SARIMAX
    except Exception as exc:
        return None, exc
    return SARIMAX, None


def _as_1d(values: float | Sequence[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values.reshape(-1)


def _poly_power(poly: np.ndarray, power: int) -> np.ndarray:
    result = np.array([1.0])
    for _ in range(int(power)):
        result = np.polynomial.polynomial.polymul(result, poly)
    return result


def _seasonal_polynomial(coefficients: np.ndarray, period: int, sign: float) -> np.ndarray:
    if period <= 1 or coefficients.size == 0:
        return np.array([1.0])
    poly = np.zeros(period * coefficients.size + 1)
    poly[0] = 1.0
    poly[period * np.arange(1, coefficients.size + 1)] = sign * coefficients
    return poly


def _validate_order(name: str, value: int) -> int:
    value = int(value)
    if value < 0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _validate_coefficients(
    name: str,
    values: Sequence[float] | np.ndarray | None,
    order: int,
) -> np.ndarray | None:
    if values is None:
        return None
    values = _as_1d(values)
    if values.size != order:
        raise ValueError(f"dimension of {name} does not match model order")
    return values


def stationary_ar(p: int, *, rng: SeedLike = None) -> np.ndarray:
    """Draw AR coefficients from the stationary region.

    The sampling scheme follows the R implementation: for orders above two it
    draws inverse roots inside the unit circle and converts them back to AR
    polynomial coefficients.
    """
    p = _validate_order("p", p)
    if p < 1:
        raise ValueError("p must be a positive integer")

    generator = get_rng(rng)
    if p == 1:
        return generator.uniform(-1.0, 1.0, size=1)
    if p == 2:
        phi2 = float(generator.uniform(-1.0, 1.0))
        phi1 = float(generator.uniform(phi2 - 1.0, 1.0 - phi2))
        return np.array([phi1, phi2])

    n_real_roots = p % 2
    inv_roots: list[complex] = [complex(x) for x in generator.uniform(-1.0, 1.0, size=n_real_roots)]
    n_complex_roots = (p - n_real_roots) // 2
    radii = generator.uniform(-1.0, 1.0, size=n_complex_roots)
    angles = generator.uniform(-np.pi, np.pi, size=n_complex_roots)
    for radius, angle in zip(radii, angles):
        root = radius * np.exp(1j * angle)
        inv_roots.extend([root, np.conjugate(root)])

    roots = np.array([1.0 / root for root in inv_roots], dtype=complex)
    poly = np.poly(roots)[::-1]
    poly = poly / poly[0]
    phi = -poly[1:].real
    return phi.astype(float)


def pi_coefficients(
    ar: float | Sequence[float] = 0.0,
    d: int = 0,
    ma: float | Sequence[float] = 0.0,
    sar: float | Sequence[float] = 0.0,
    D: int = 0,
    sma: float | Sequence[float] = 0.0,
    m: int = 1,
    tol: float = 1e-7,
    *,
    n_terms: int = 500,
) -> np.ndarray:
    """Convert SARIMA coefficients to AR recursion coefficients.

    This mirrors the exported R helper and is useful for the MAR simulation
    internals. Polynomial coefficients are stored in ascending lag order.
    """
    ar = _as_1d(ar)
    ma = _as_1d(ma)
    sar = _as_1d(sar)
    sma = _as_1d(sma)
    m = int(m)

    ar_poly = np.polynomial.polynomial.polymul(
        np.concatenate(([1.0], -ar)),
        _poly_power(np.array([1.0, -1.0]), int(d)),
    )

    if m > 1:
        seasonal_ar = _seasonal_polynomial(sar, m, sign=-1.0)
        seasonal_diff = np.zeros(m + 1)
        seasonal_diff[0] = 1.0
        seasonal_diff[m] = -1.0
        sar_poly = np.polynomial.polynomial.polymul(
            seasonal_ar,
            _poly_power(seasonal_diff, int(D)),
        )
    else:
        sar_poly = np.array([1.0])

    ma_poly = np.concatenate(([1.0], ma))
    if m > 1:
        sma_poly = _seasonal_polynomial(sma, m, sign=1.0)
    else:
        sma_poly = np.array([1.0])

    theta = -np.polynomial.polynomial.polymul(ma_poly, sma_poly)[1:]
    if theta.size == 0:
        theta = np.array([0.0])

    phi = -np.concatenate((
        np.polynomial.polynomial.polymul(ar_poly, sar_poly)[1:],
        np.zeros(int(n_terms)),
    ))
    q = len(theta)
    pie = np.zeros(q + int(n_terms) + 1)
    pie[q] = 1.0
    for j in range(int(n_terms)):
        lagged = pie[j + 1:j + q + 1][::-1]
        pie[j + q + 1] = -phi[j] - float(np.dot(theta, lagged))

    coeffs = -pie[q + 1:q + int(n_terms) + 1]
    nz = np.flatnonzero(np.abs(coeffs) > tol)
    return coeffs[: nz[-1] + 1] if nz.size else np.array([], dtype=float)


@dataclass(frozen=True)
class ARIMAModel:
    """Gaussian SARIMA model.

    Parameters follow the R package's `arima_model()` function. `sigma` is the
    innovation standard deviation.
    """

    frequency: int = 1
    order: tuple[int, int, int] = (0, 0, 0)
    seasonal_order: tuple[int, int, int] = (0, 0, 0)
    constant: float = 0.0
    phi: np.ndarray | None = None
    theta: np.ndarray | None = None
    Phi: np.ndarray | None = None
    Theta: np.ndarray | None = None
    sigma: float = 1.0
    backend: str = "auto"

    def __post_init__(self) -> None:
        object.__setattr__(self, "phi", np.asarray([] if self.phi is None else self.phi, dtype=float))
        object.__setattr__(self, "theta", np.asarray([] if self.theta is None else self.theta, dtype=float))
        object.__setattr__(self, "Phi", np.asarray([] if self.Phi is None else self.Phi, dtype=float))
        object.__setattr__(self, "Theta", np.asarray([] if self.Theta is None else self.Theta, dtype=float))
        backend = self.backend.lower()
        if backend not in {"auto", "statsmodels", "numpy"}:
            raise ValueError("backend must be one of 'auto', 'statsmodels', or 'numpy'")
        object.__setattr__(self, "backend", backend)

    @property
    def p(self) -> int:
        return self.order[0]

    @property
    def d(self) -> int:
        return self.order[1]

    @property
    def q(self) -> int:
        return self.order[2]

    @property
    def P(self) -> int:
        return self.seasonal_order[0]

    @property
    def D(self) -> int:
        return self.seasonal_order[1]

    @property
    def Q(self) -> int:
        return self.seasonal_order[2]

    def ar_polynomial(self, *, include_differences: bool = True) -> np.ndarray:
        """Return the AR-side polynomial in ascending lag order."""
        nonseasonal = np.concatenate(([1.0], -self.phi))
        seasonal = _seasonal_polynomial(self.Phi, self.frequency, sign=-1.0)
        ar_poly = np.polynomial.polynomial.polymul(nonseasonal, seasonal)
        if include_differences and self.d:
            ar_poly = np.polynomial.polynomial.polymul(
                ar_poly,
                _poly_power(np.array([1.0, -1.0]), self.d),
            )
        if include_differences and self.D and self.frequency > 1:
            seasonal_diff = np.zeros(self.frequency + 1)
            seasonal_diff[0] = 1.0
            seasonal_diff[self.frequency] = -1.0
            ar_poly = np.polynomial.polynomial.polymul(
                ar_poly,
                _poly_power(seasonal_diff, self.D),
            )
        return ar_poly

    def ma_polynomial(self) -> np.ndarray:
        """Return the MA-side polynomial in ascending lag order."""
        nonseasonal = np.concatenate(([1.0], self.theta))
        seasonal = _seasonal_polynomial(self.Theta, self.frequency, sign=1.0)
        return np.polynomial.polynomial.polymul(nonseasonal, seasonal)

    def simulate(
        self,
        n: int = 100,
        *,
        n_start: int = 100,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Simulate one sample path from the SARIMA model."""
        n = _validate_order("n", n)
        n_start = _validate_order("n_start", n_start)
        if n == 0:
            return np.array([], dtype=float)
        if self.sigma <= 0:
            raise ValueError("sigma must be positive")

        if self.backend in {"auto", "statsmodels"}:
            model_class, error = _load_statsmodels_sarimax()
            if model_class is not None:
                return self._simulate_statsmodels(
                    model_class,
                    n=n,
                    n_start=n_start,
                    rng=rng,
                )
            if self.backend == "statsmodels":
                raise ImportError(
                    "statsmodels could not be imported for ARIMA simulation. "
                    "Install the Python statsmodels extra with "
                    "`python -m pip install -e './python[statsmodels]'`."
                ) from error

        return self._simulate_numpy(n=n, n_start=n_start, rng=rng)

    def _simulate_statsmodels(
        self,
        model_class: type[Any],
        *,
        n: int,
        n_start: int = 100,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Simulate with statsmodels' SARIMAX backend."""
        generator = get_rng(rng)
        seasonal_period = self.frequency if any(self.seasonal_order) else 0
        trend = "c" if self.constant != 0 else "n"
        endog_length = max(10, self.frequency * 2)
        model = model_class(
            np.zeros(endog_length, dtype=float),
            order=self.order,
            seasonal_order=(*self.seasonal_order, seasonal_period),
            trend=trend,
            enforce_stationarity=False,
            enforce_invertibility=False,
        )

        params: list[float] = []
        if self.constant != 0:
            params.append(float(self.constant))
        params.extend(np.asarray(self.phi, dtype=float).tolist())
        params.extend(np.asarray(self.theta, dtype=float).tolist())
        params.extend(np.asarray(self.Phi, dtype=float).tolist())
        params.extend(np.asarray(self.Theta, dtype=float).tolist())
        params.append(float(self.sigma) ** 2)

        total = n + n_start
        simulated = model.simulate(
            np.asarray(params, dtype=float),
            nsimulations=total,
            initial_state=np.zeros(model.k_states, dtype=float),
            random_state=generator,
        )
        return np.asarray(simulated, dtype=float).reshape(-1)[-n:]

    def _simulate_numpy(
        self,
        n: int = 100,
        *,
        n_start: int = 100,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Fallback SARIMA simulation used when statsmodels is unavailable."""
        generator = get_rng(rng)
        ar_poly = self.ar_polynomial(include_differences=True)
        ma_poly = self.ma_polynomial()
        ar_lag = ar_poly.size - 1
        ma_lag = ma_poly.size - 1
        start = max(ar_lag, ma_lag)
        total = n + n_start + start
        innovations = generator.normal(loc=0.0, scale=self.sigma, size=total)
        y = np.zeros(total, dtype=float)

        for t in range(start, total):
            ar_part = 0.0
            if ar_lag:
                ar_part = -float(np.dot(ar_poly[1:], y[t - np.arange(1, ar_lag + 1)]))
            ma_part = 0.0
            if ma_lag:
                ma_part = float(np.dot(ma_poly[1:], innovations[t - np.arange(1, ma_lag + 1)]))
            y[t] = self.constant + ar_part + innovations[t] + ma_part

        return y[-n:]


def arima_model(
    frequency: int = 1,
    p: int | None = None,
    d: int | None = None,
    q: int | None = None,
    P: int | None = None,
    D: int | None = None,
    Q: int | None = None,
    constant: float | None = None,
    phi: Sequence[float] | np.ndarray | None = None,
    theta: Sequence[float] | np.ndarray | None = None,
    Phi: Sequence[float] | np.ndarray | None = None,
    Theta: Sequence[float] | np.ndarray | None = None,
    sigma: float | None = None,
    backend: str = "auto",
    *,
    rng: SeedLike = None,
) -> ARIMAModel:
    """Specify a Gaussian SARIMA model and return a simulatable object."""
    generator = get_rng(rng)
    frequency = int(frequency)
    if frequency < 1:
        raise ValueError("frequency must be a positive integer")

    if sigma is None:
        sigma = float(generator.uniform(1.0, 5.0))
    if sigma <= 0:
        raise ValueError("sigma must be positive")

    if p is None:
        p = len(_as_1d(phi)) if phi is not None else int(generator.choice([0, 1, 2, 3]))
    if d is None:
        d = int(generator.choice([0, 1, 2]))
    if q is None:
        q = len(_as_1d(theta)) if theta is not None else int(generator.choice([0, 1, 2, 3]))
    p = _validate_order("p", p)
    d = _validate_order("d", d)
    q = _validate_order("q", q)

    if frequency > 1:
        if P is None:
            P = len(_as_1d(Phi)) if Phi is not None else int(generator.choice([0, 1, 2]))
        if D is None:
            D = 0 if d == 2 else int(generator.choice([0, 1]))
        if Q is None:
            Q = len(_as_1d(Theta)) if Theta is not None else int(generator.choice([0, 1, 2]))
    else:
        seasonal_orders = [value for value in (P, D, Q) if value is not None]
        supplied_seasonal_coefficients = (
            Phi is not None and _as_1d(Phi).size > 0
        ) or (
            Theta is not None and _as_1d(Theta).size > 0
        )
        if any(int(value) != 0 for value in seasonal_orders) or supplied_seasonal_coefficients:
            raise ValueError("seasonal arguments require frequency greater than 1")
        P = D = Q = 0

    P = _validate_order("P", P)
    D = _validate_order("D", D)
    Q = _validate_order("Q", Q)
    if d + D > 2:
        raise ValueError("d + D must be no greater than 2")

    phi_values = _validate_coefficients("phi", phi, p)
    theta_values = _validate_coefficients("theta", theta, q)
    Phi_values = _validate_coefficients("Phi", Phi, P)
    Theta_values = _validate_coefficients("Theta", Theta, Q)

    if phi_values is None:
        phi_values = stationary_ar(p, rng=generator) if p > 0 else np.array([], dtype=float)
    if theta_values is None:
        theta_values = -stationary_ar(q, rng=generator) if q > 0 else np.array([], dtype=float)
    if Phi_values is None:
        Phi_values = stationary_ar(P, rng=generator) if P > 0 else np.array([], dtype=float)
    if Theta_values is None:
        Theta_values = -stationary_ar(Q, rng=generator) if Q > 0 else np.array([], dtype=float)

    if constant is None:
        constant = float(generator.uniform(-3.0, 3.0)) if d + D <= 1 else 0.0

    return ARIMAModel(
        frequency=frequency,
        order=(p, d, q),
        seasonal_order=(P, D, Q),
        constant=float(constant),
        phi=phi_values,
        theta=theta_values,
        Phi=Phi_values,
        Theta=Theta_values,
        sigma=float(sigma),
        backend=backend,
    )
