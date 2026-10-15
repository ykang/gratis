"""Mixture autoregressive model structure and simulation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ._random import SeedLike, get_rng
from .arima import _validate_order, pi_coefficients, stationary_ar
from .mixture import _normalise_weights


def _as_int_vector(name: str, values: int | Sequence[int], k: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 1:
        arr = np.repeat(arr, k)
    if arr.size != k:
        raise ValueError(f"dimension of {name} does not match other components")
    if np.any(~np.isfinite(arr)) or np.any(arr < 0) or np.any(arr != np.floor(arr)):
        raise ValueError(f"{name} must contain non-negative integers")
    return arr.astype(int)


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
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


@dataclass(frozen=True)
class MARModel:
    """Store a fully expanded mixture autoregressive model.

    Use ``mar_model`` to validate inputs and construct consistent metadata and
    recursion coefficients. This low-level dataclass copies supplied arrays and
    makes them read-only; it does not derive ``ar`` from the other fields.

    Attributes
    ----------
    k : int
        Number of mixture components.
    p, d, P, D : numpy.ndarray, shape (k,)
        Non-seasonal AR, ordinary difference, seasonal AR and seasonal
        difference orders, respectively.
    phi, Phi : numpy.ndarray
        Non-seasonal and seasonal coefficient matrices, with components in
        columns and zero padding for lower-order components.
    constants : numpy.ndarray, shape (k,)
        Component recursion intercepts, not unconditional means.
    ar : numpy.ndarray, shape (max_order + 1, k)
        Expanded coefficients. Row zero contains the intercept; subsequent
        rows contain the coefficients of lags one through ``max_order``.
    sigmas : numpy.ndarray, shape (k,)
        Positive innovation variances, not standard deviations.
    weights : numpy.ndarray, shape (k,)
        Normalized component probabilities.
    m : numpy.ndarray, shape (k,)
        Positive integer seasonal period for each component.

    Notes
    -----
    Construct a new model with ``mar_model`` to change parameters. Stationary
    component AR polynomials alone do not ensure mixture stability, and
    integrated components are permitted.
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

    def __post_init__(self) -> None:
        for name in ("p", "d", "phi", "P", "D", "Phi", "constants", "ar",
                     "sigmas", "weights", "m"):
            values = np.array(getattr(self, name), copy=True)
            values.setflags(write=False)
            object.__setattr__(self, name, values)

    @property
    def n_components(self) -> int:
        """Number of mixture components."""
        return int(self.k)

    @property
    def seasonal_periods(self) -> np.ndarray:
        """Read-only seasonal-period vector, with one entry per component."""
        return self.m

    def simulate(
        self,
        n: int = 100,
        *,
        n_start: int = 100,
        rng: SeedLike = None,
    ) -> np.ndarray:
        """Simulate from zero states, discarding n_start observations.

        The default 100-step burn-in may be insufficient near a unit root.
        Stationary component polynomials do not guarantee mixture stability.

        Parameters
        ----------
        n : int, default 100
            Non-negative number of returned observations.
        n_start : int, default 100
            Non-negative number of recursively generated observations to discard
            after zero prehistory. Zero disables burn-in.
        rng : int or numpy.random.Generator, optional
            Seed or generator. Existing generators advance in place; ``None``
            creates a new independently seeded generator.

        Returns
        -------
        numpy.ndarray
            Path of shape ``(n,)``, without calendar metadata or standardization.
            A component is drawn independently at each step, and its conditional
            mean uses the common observed past, not a private component path.

        Raises
        ------
        ValueError
            If ``n`` or ``n_start`` is not a non-negative integer.

        See Also
        --------
        gratis.generate : Generate several paths with successive random draws.
        """
        n = _validate_order("n", n)
        n_start = _validate_order("n_start", n_start)
        if n == 0:
            return np.array([], dtype=float)

        generator = get_rng(rng)
        order = self.ar.shape[0] - 1
        total = n + n_start + order
        components = generator.choice(self.k, size=total, replace=True, p=self.weights)
        selected_variances = self.sigmas[components]
        innovations = generator.normal(loc=0.0, scale=np.sqrt(selected_variances), size=total)
        if order == 0:
            return (self.ar[0, components] + innovations)[n_start:]
        if self.k == 1:
            from statsmodels.tsa.arima_process import arma_generate_sample

            noise = self.ar[0, 0] + innovations[order:]
            y = arma_generate_sample(np.r_[1., -self.ar[1:, 0]], [1.],
                                     nsample=noise.size, distrvs=lambda size: noise)
            return np.asarray(y)[n_start:]
        y = np.zeros(total, dtype=float)
        lags = [np.flatnonzero(self.ar[1:, j]) + 1 for j in range(self.k)]
        coefficients = [self.ar[lags[j], j] for j in range(self.k)]

        for t in range(order, total):
            component = components[t]
            y[t] = (self.ar[0, component]
                    + float(np.dot(coefficients[component], y[t - lags[component]]))
                    + innovations[t])

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
    """Construct a mixture of Gaussian autoregressive components.

    Parameters
    ----------
    k : int, optional
        Positive component count. Inferred from supplied arrays or seasonal
        periods where possible; otherwise sampled from ``{1, 2, 3, 4, 5}``.
    p, P : int or array_like, optional
        Non-negative non-seasonal and seasonal AR orders. Scalars apply to all
        components. Supplied ``phi`` or ``Phi`` determines the corresponding
        orders instead. Missing orders are sampled from ``{0, 1, 2, 3}`` and
        ``{0, 1, 2}``, respectively.
    d, D : int or array_like, optional
        Non-negative ordinary and seasonal difference orders. Scalars apply to
        all components; missing values are sampled from ``{0, 1, 2}`` and
        ``{0, 1}``. Seasonal differences are zero when ``d=2`` or the period
        is one. All seasonal terms are ignored when every period is one.
    phi, Phi : array_like, optional
        Finite coefficient matrices with one column per component. A vector
        specifies one component; pad with zeros for different component orders.
        Missing coefficients are sampled from the stationary AR region.
    constants : array_like, shape (k,), optional
        Recursion intercepts, sampled uniformly from ``(-3, 3)`` when omitted.
        Forced to zero for components with ``d + D > 1``. These are not ARIMA
        mean/drift parameters and are not replicated from a scalar for ``k>1``.
    sigmas : array_like, shape (k,), optional
        Positive finite innovation variances, not standard deviations. Sampled
        uniformly from ``(1, 5)`` when omitted. Python does not accept the R
        implementation's GARCH specifications.
    weights : array_like, shape (k,), optional
        Finite, non-negative weights with at least one positive value. Supplied
        weights are normalized; missing weights are independent uniform draws
        on ``(0, 1)`` before normalization.
    seasonal_periods : int or array_like, default 1
        Positive integer period per component. A scalar is repeated ``k`` times.
    rng : int or numpy.random.Generator, optional
        Seed or generator used to sample missing parameters. An existing
        generator advances in place; ``None`` creates a new one. A constructor
        seed is not stored as the seed for future simulations.

    Returns
    -------
    MARModel
        Model with copied, read-only arrays and expanded recursion coefficients.

    Raises
    ------
    ValueError
        If dimensions, orders, periods, weights or numerical parameters are
        invalid. Explicit coefficients need not define a stationary process.

    See Also
    --------
    MARModel.simulate : Simulate from zero states with configurable burn-in.
    pi_coefficients : Expand seasonal and differencing polynomials.

    Examples
    --------
    >>> model = mar_model(phi=[[0.4, -0.2]], d=0, constants=[0.5, -0.2],
    ...                   sigmas=[1, 2], weights=[0.7, 0.3])
    >>> model.simulate(12, rng=1).shape
    (12,)
    """
    generator = get_rng(rng)

    seasonal_periods_arr = np.asarray(seasonal_periods, dtype=float).reshape(-1)
    if (seasonal_periods_arr.size == 0 or np.any(~np.isfinite(seasonal_periods_arr))
            or np.any(seasonal_periods_arr < 1)
            or np.any(seasonal_periods_arr != np.floor(seasonal_periods_arr))):
        raise ValueError("seasonal_periods must contain positive integers")
    seasonal_periods_arr = seasonal_periods_arr.astype(int)

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
    k = _validate_order("k", k)
    if k < 1:
        raise ValueError("k must be positive")

    if weights is None:
        weights_arr = generator.uniform(0.0, 1.0, size=k)
    else:
        weights_arr = np.asarray(weights, dtype=float).reshape(-1)
        if weights_arr.size != k:
            raise ValueError("dimension of weights does not match other components")
    weights_arr = _normalise_weights(weights_arr, k)

    if constants_arr is None:
        constants_arr = generator.uniform(-3.0, 3.0, size=k)
    elif constants_arr.size != k:
        raise ValueError("dimension of constants does not match other components")
    if np.any(~np.isfinite(constants_arr)):
        raise ValueError("constants must be finite")

    if sigmas_arr is None:
        sigmas_arr = generator.uniform(1.0, 5.0, size=k)
    elif sigmas_arr.size != k:
        raise ValueError("dimension of sigmas does not match other components")
    if np.any(~np.isfinite(sigmas_arr)) or np.any(sigmas_arr <= 0):
        raise ValueError("sigmas must be positive and finite")

    # Validate supplied values even when a period or explicit coefficient array
    # makes an order redundant in the expanded model.
    for name, values in (("p", p), ("d", d), ("P", P), ("D", D)):
        if values is not None:
            _as_int_vector(name, values, k)
    if Phi_matrix is not None:
        Phi_matrix = _as_coefficient_matrix("Phi", Phi_matrix, k)

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
