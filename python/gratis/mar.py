"""Mixture autoregressive model structure and simulation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ._random import SeedLike, get_rng
from .arima import pi_coefficients, stationary_ar
from .mixture import _normalise_weights


def _as_int_vector(name: str, values: int | Sequence[int], k: int) -> np.ndarray:
    arr = np.asarray(values, dtype=int).reshape(-1)
    if arr.size == 1:
        arr = np.repeat(arr, k)
    if arr.size != k:
        raise ValueError(f"dimension of {name} does not match other components")
    if np.any(arr < 0):
        raise ValueError(f"{name} must be non-negative")
    return arr


def _as_float_vector(name: str, values: float | Sequence[float], k: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 1 and k > 1:
        arr = np.repeat(arr, k)
    if arr.size != k:
        raise ValueError(f"dimension of {name} does not match other components")
    return arr


def _as_coefficient_matrix(
    name: str,
    values: Sequence[float] | np.ndarray,
    k: int,
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2 or arr.shape[1] != k:
        raise ValueError(f"dimension of {name} does not match other components")
    return arr


@dataclass(frozen=True)
class MARModel:
    """Mixture autoregressive model.

    The expanded `ar` matrix has one column per mixture component. Row zero is
    the component constant and later rows are lag coefficients.
    """

    k: int
    p: np.ndarray
    d: np.ndarray
    phi: np.ndarray
    P: np.ndarray
    D: np.ndarray
    Phi: np.ndarray
    constants: np.ndarray
    ar: np.ndarray
    sigmas: np.ndarray
    weights: np.ndarray
    m: np.ndarray

    @property
    def n_components(self) -> int:
        return int(self.k)

    @property
    def seasonal_periods(self) -> np.ndarray:
        return self.m

    def simulate(
        self,
        n: int = 100,
        *,
        n_start: int = 100,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Simulate one sample path from the MAR model."""
        n = int(n)
        n_start = int(n_start)
        if n < 0:
            raise ValueError("n must be non-negative")
        if n_start < 0:
            raise ValueError("n_start must be non-negative")
        if n == 0:
            return np.array([], dtype=float)

        generator = get_rng(rng)
        order = self.ar.shape[0] - 1
        total = n + n_start + order
        components = generator.choice(self.k, size=total, replace=True, p=self.weights)
        selected_variances = self.sigmas[components]
        innovations = generator.normal(loc=0.0, scale=np.sqrt(selected_variances), size=total)
        y = np.zeros(total, dtype=float)

        for t in range(order, total):
            component = components[t]
            coeffs = self.ar[:, component]
            lagged = y[t - np.arange(1, order + 1)] if order else np.array([], dtype=float)
            y[t] = coeffs[0] + float(np.dot(coeffs[1:], lagged)) + innovations[t]

        return y[n_start + order :]

    def __str__(self) -> str:
        lines = [f"Mixture AR model with {self.k} components:"]
        for index in range(self.k):
            model = f"ARIMA({self.p[index]},{self.d[index]},0)"
            if self.P[index] > 0 or self.D[index] > 0:
                model += f"({self.P[index]},{self.D[index]},0)[{self.m[index]}]"
            lines.append(f"    {model} with weight {self.weights[index]:.2f}")
        return "\n".join(lines)


def mar_model(
    k: int | None = None,
    p: int | Sequence[int] | None = None,
    d: int | Sequence[int] | None = None,
    phi: Sequence[float] | np.ndarray | None = None,
    P: int | Sequence[int] | None = None,
    D: int | Sequence[int] | None = None,
    Phi: Sequence[float] | np.ndarray | None = None,
    constants: float | Sequence[float] | None = None,
    sigmas: float | Sequence[float] | None = None,
    weights: Sequence[float] | None = None,
    seasonal_periods: int | Sequence[int] = 1,
    *,
    rng: SeedLike = None,
) -> MARModel:
    """Specify a mixture autoregressive model.

    This is the Python analogue of the R `mar_model()` constructor. Numeric
    `sigmas` are treated as innovation variances, matching the R simulation
    method.
    """
    generator = get_rng(rng)

    seasonal_periods_arr = np.asarray(seasonal_periods, dtype=int).reshape(-1)
    if seasonal_periods_arr.size == 0 or np.any(seasonal_periods_arr < 1):
        raise ValueError("seasonal_periods must contain positive integers")

    phi_matrix = None if phi is None else np.asarray(phi, dtype=float)
    Phi_matrix = None if Phi is None else np.asarray(Phi, dtype=float)
    sigmas_arr = None if sigmas is None else np.asarray(sigmas, dtype=float).reshape(-1)
    constants_arr = None if constants is None else np.asarray(constants, dtype=float).reshape(-1)

    if k is None:
        if weights is not None:
            k = len(np.asarray(weights).reshape(-1))
        elif phi_matrix is not None:
            k = phi_matrix.shape[1] if phi_matrix.ndim == 2 else 1
        elif Phi_matrix is not None:
            k = Phi_matrix.shape[1] if Phi_matrix.ndim == 2 else 1
        elif sigmas_arr is not None:
            k = sigmas_arr.size
        elif constants_arr is not None:
            k = constants_arr.size
        elif seasonal_periods_arr.size == 1:
            k = int(generator.integers(1, 6))
        else:
            k = seasonal_periods_arr.size
    k = int(k)
    if k < 1:
        raise ValueError("k must be positive")

    if weights is None:
        weights_arr = generator.uniform(0.0, 1.0, size=k)
    else:
        weights_arr = np.asarray(weights, dtype=float).reshape(-1)
        if weights_arr.size != k:
            raise ValueError("dimension of weights does not match other components")
        if np.any(weights_arr <= 0):
            raise ValueError("weights must be positive")
    weights_arr = _normalise_weights(weights_arr, k)

    if constants_arr is None:
        constants_arr = generator.uniform(-3.0, 3.0, size=k)
    elif constants_arr.size != k:
        raise ValueError("dimension of constants does not match other components")

    if sigmas_arr is None:
        sigmas_arr = generator.uniform(1.0, 5.0, size=k)
    elif sigmas_arr.size != k:
        raise ValueError("dimension of sigmas does not match other components")
    if np.any(sigmas_arr <= 0):
        raise ValueError("sigmas must be positive")

    if seasonal_periods_arr.size == 1:
        seasonal_periods_arr = np.repeat(seasonal_periods_arr, k)
    elif seasonal_periods_arr.size != k:
        raise ValueError("dimension of seasonal_periods does not match other components")

    if p is None and phi_matrix is None:
        p_arr = generator.choice([0, 1, 2, 3], size=k, replace=True)
    elif phi_matrix is not None:
        phi_matrix = _as_coefficient_matrix("phi", phi_matrix, k)
        p_arr = np.repeat(phi_matrix.shape[0], k)
    else:
        p_arr = _as_int_vector("p", p, k)

    if phi_matrix is None:
        max_p = int(p_arr.max(initial=0))
        phi_matrix = np.zeros((max_p, k), dtype=float)
        for index, order in enumerate(p_arr):
            if order > 0:
                phi_matrix[:order, index] = stationary_ar(int(order), rng=generator)

    if d is None:
        d_arr = generator.choice([0, 1, 2], size=k, replace=True)
    else:
        d_arr = _as_int_vector("d", d, k)

    if np.any(seasonal_periods_arr > 1):
        if P is None and Phi_matrix is None:
            P_arr = generator.choice([0, 1, 2], size=k, replace=True)
        elif Phi_matrix is not None:
            Phi_matrix = _as_coefficient_matrix("Phi", Phi_matrix, k)
            P_arr = np.repeat(Phi_matrix.shape[0], k)
        else:
            P_arr = _as_int_vector("P", P, k)

        if Phi_matrix is None:
            P_arr = np.asarray(P_arr, dtype=int)
            P_arr[seasonal_periods_arr == 1] = 0
            max_P = int(P_arr.max(initial=0))
            Phi_matrix = np.zeros((max_P, k), dtype=float)
            for index, order in enumerate(P_arr):
                if order > 0:
                    Phi_matrix[:order, index] = stationary_ar(int(order), rng=generator)

        if D is None:
            D_arr = generator.choice([0, 1], size=k, replace=True)
        else:
            D_arr = _as_int_vector("D", D, k)
        D_arr = np.asarray(D_arr, dtype=int)
        D_arr[(seasonal_periods_arr == 1) | (d_arr == 2)] = 0
    else:
        P_arr = np.zeros(k, dtype=int)
        D_arr = np.zeros(k, dtype=int)
        Phi_matrix = np.zeros((0, k), dtype=float)

    constants_arr = constants_arr.astype(float, copy=True)
    constants_arr[d_arr + D_arr > 1] = 0.0

    max_order = int(np.max(d_arr + seasonal_periods_arr * D_arr + p_arr + seasonal_periods_arr * P_arr))
    ar_matrix = np.zeros((max_order + 1, k), dtype=float)
    for index in range(k):
        pi_weights = np.concatenate(
            (
                np.array([constants_arr[index]], dtype=float),
                pi_coefficients(
                    ar=phi_matrix[:, index] if phi_matrix.size else np.array([], dtype=float),
                    sar=Phi_matrix[:, index] if Phi_matrix.size else np.array([], dtype=float),
                    d=int(d_arr[index]),
                    D=int(D_arr[index]),
                    m=int(seasonal_periods_arr[index]),
                ),
            )
        )
        ar_matrix[: pi_weights.size, index] = pi_weights

    nonzero = np.flatnonzero(np.sum(np.abs(ar_matrix), axis=1) > 0)
    if nonzero.size:
        ar_matrix = ar_matrix[: nonzero[-1] + 1, :]
    else:
        ar_matrix = ar_matrix[:1, :]

    return MARModel(
        k=k,
        p=np.asarray(p_arr, dtype=int),
        d=np.asarray(d_arr, dtype=int),
        phi=np.asarray(phi_matrix, dtype=float),
        P=np.asarray(P_arr, dtype=int),
        D=np.asarray(D_arr, dtype=int),
        Phi=np.asarray(Phi_matrix, dtype=float),
        constants=constants_arr,
        ar=ar_matrix,
        sigmas=sigmas_arr.astype(float),
        weights=weights_arr,
        m=seasonal_periods_arr.astype(int),
    )
