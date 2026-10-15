"""ARIMA utilities and simulation models for the Python gratis implementation."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import warnings

import numpy as np

from ._random import SeedLike, get_rng


def _load_statsmodels_arma() -> tuple[Callable | None, Exception | None]:
    try:
        from statsmodels.tsa.arima_process import arma_generate_sample
    except Exception as exc:
        return None, exc
    return arma_generate_sample, None


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
    if np.ndim(value) != 0 or not np.isfinite(value) or value != int(value):
        raise ValueError(f"{name} must be an integer")
    value = int(value)
    if value < 0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _validate_sigma(sigma: float) -> None:
    if np.ndim(sigma) != 0 or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be positive and finite")
    limits = np.finfo(float)
    if sigma < np.sqrt(limits.tiny) or sigma > np.sqrt(limits.max):
        raise ValueError("sigma squared must be a positive finite normal variance")


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
    if np.any(~np.isfinite(values)):
        raise ValueError(f"{name} must be finite")
    return values


def stationary_ar(p: int, *, rng: SeedLike = None) -> np.ndarray:
    """Draw AR coefficients from the stationary region.

    The sampling scheme follows the R implementation: for orders above two it
    draws inverse roots inside the unit circle, allowing every possible number
    of real roots, and converts them back to AR polynomial coefficients. This
    has full support but is not uniform over the stationary region.

    Parameters
    ----------
    p : int
        Positive autoregressive order.
    rng : int or numpy.random.Generator, optional
        Seed or generator. An existing generator advances in place; ``None``
        creates a new independently seeded generator.

    Returns
    -------
    numpy.ndarray
        Shape ``(p,)``, using the convention
        ``y[t] = phi[0]*y[t-1] + ... + phi[p-1]*y[t-p] + error[t]``.

    Raises
    ------
    ValueError
        If ``p`` is not a positive integer.

    Examples
    --------
    >>> stationary_ar(3, rng=1).shape
    (3,)
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

    n_real_roots = int(generator.choice(np.arange(p % 2, p + 1, 2)))
    inv_roots: list[complex] = [complex(x) for x in generator.uniform(-1.0, 1.0, size=n_real_roots)]
    n_complex_roots = (p - n_real_roots) // 2
    radii = generator.uniform(-1.0, 1.0, size=n_complex_roots)
    angles = generator.uniform(-np.pi, np.pi, size=n_complex_roots)
    for radius, angle in zip(radii, angles):
        root = radius * np.exp(1j * angle)
        inv_roots.extend([root, np.conjugate(root)])

    # Characteristic polynomial in inverse roots avoids division by tiny roots.
    poly = np.poly(inv_roots)
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

    Polynomial coefficients are stored in ascending lag order, with positive
    MA signs as in R's ARIMA models. Without MA terms the complete finite AR
    polynomial is returned. With MA terms, n_terms is increased to include all
    polynomial degrees and a warning signals an unconverged approximation.

    Parameters
    ----------
    ar, ma : float or array_like, default 0
        Finite non-seasonal AR and MA coefficients, excluding leading ones.
    sar, sma : float or array_like, default 0
        Seasonal AR and MA coefficients. Ignored when ``m=1``.
    d, D : int, default 0
        Non-negative ordinary and seasonal difference orders. Seasonal
        differencing is ignored when ``m=1``.
    m : int, default 1
        Positive seasonal period; one denotes non-seasonal data.
    tol : float, default 1e-7
        Non-negative truncation threshold. Coefficients through the last lag
        whose absolute value exceeds this threshold are retained.
    n_terms : int, default 500
        Positive lag limit for an infinite AR approximation with MA terms.
        Increased as needed to include the full polynomial degrees.

    Returns
    -------
    numpy.ndarray
        AR recursion coefficients in increasing lag order, without an intercept
        or leading one. Complete cancellation returns an empty array.

    Warns
    -----
    RuntimeWarning
        If the MA approximation has not decayed below ``tol``. Increase
        ``n_terms``; noninvertible MA terms may not yield a convergent expansion.

    Examples
    --------
    >>> pi_coefficients(ar=[0.5]).tolist()
    [0.5]
    """
    ar = _as_1d(ar)
    ma = _as_1d(ma)
    sar = _as_1d(sar)
    sma = _as_1d(sma)
    if any(np.any(~np.isfinite(values)) for values in (ar, ma, sar, sma)):
        raise ValueError("AR and MA coefficients must be finite")
    m = _validate_order("m", m)
    d = _validate_order("d", d)
    D = _validate_order("D", D)
    n_terms = _validate_order("n_terms", n_terms)
    if m < 1 or n_terms < 1:
        raise ValueError("m and n_terms must be positive integers")
    if not np.isfinite(tol) or tol < 0:
        raise ValueError("tol must be non-negative and finite")

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

    numerator = np.polynomial.polynomial.polymul(ar_poly, sar_poly)
    denominator = np.polynomial.polynomial.polymul(ma_poly, sma_poly)
    q = denominator.size - 1
    if q == 0:
        coeffs = -numerator[1:]
    else:
        from statsmodels.tsa.arima_process import arma2ar

        n = max(n_terms, numerator.size - 1, q)
        coeffs = -arma2ar(numerator, denominator, lags=n + 1)[1:]
        if np.any(~np.isfinite(coeffs)) or np.any(np.abs(coeffs[-q:]) > tol):
            warnings.warn(
                "MA approximation has not decayed below tol; increase n_terms",
                RuntimeWarning,
                stacklevel=2,
            )
    nz = np.flatnonzero(np.abs(coeffs) > tol)
    return coeffs[: nz[-1] + 1] if nz.size else np.array([], dtype=float)


@dataclass(frozen=True)
class ARIMAModel:
    """Store a fully specified Gaussian SARIMA model.

    Parameters
    ----------
    frequency : int, default 1
        Positive seasonal period. Seasonal terms require a value above one.
    order, seasonal_order : tuple of three ints, default (0, 0, 0)
        Orders ``(p, d, q)`` and ``(P, D, Q)``. All are non-negative and
        ``d + D <= 2``.
    constant : float, default 0
        Mean when ``d + D == 0`` or linear drift slope when ``d + D == 1``.
        Must be zero with two differences; this is not the MAR intercept.
    phi, theta, Phi, Theta : array_like, optional
        Finite non-seasonal and seasonal AR/MA coefficients of the corresponding
        orders. Missing arrays are allowed only for zero-order components.
        The MA lag polynomial has positive signs: ``1 + theta(B)``.
    sigma : float, default 1
        Innovation standard deviation. Its square must be a positive finite
        normal floating-point variance, matching the R implementation.
    backend : {"auto", "statsmodels", "numpy"}, default "auto"
        ``auto`` selects statsmodels' ARMA filter, without a silent fallback.
        ``numpy`` selects an explicit recursion for debugging.

    Notes
    -----
    Use ``arima_model`` to sample omitted parameters instead of supplying every
    coefficient. Coefficient arrays are copied at construction but remain
    writable. Stationarity of the AR part is checked before simulation.

    Automatic burn-in applies to the stationary ARMA part before integration.
    It uses the same adaptive impulse-energy estimate as R, not exact sampling
    from the stationary distribution.

    See Also
    --------
    arima_model : Construct a model with randomly selected missing parameters.
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
        frequency = _validate_order("frequency", self.frequency)
        if frequency < 1:
            raise ValueError("frequency must be a positive integer")
        object.__setattr__(self, "frequency", frequency)
        for name in ("order", "seasonal_order"):
            values = getattr(self, name)
            if np.ndim(values) != 1 or len(values) != 3:
                raise ValueError(f"{name} must contain three non-negative integers")
            object.__setattr__(self, name, tuple(_validate_order(name, value) for value in values))
        if frequency == 1 and any(self.seasonal_order):
            raise ValueError("seasonal arguments require frequency greater than 1")
        if self.d + self.D > 2:
            raise ValueError("d + D must be no greater than 2")
        for name, order in (("phi", self.p), ("theta", self.q), ("Phi", self.P), ("Theta", self.Q)):
            values = getattr(self, name)
            values = _validate_coefficients(name, [] if values is None else values, order)
            object.__setattr__(self, name, values.copy())
        for name in ("constant", "sigma"):
            value = getattr(self, name)
            if np.ndim(value) != 0 or not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
        _validate_sigma(self.sigma)
        if self.d + self.D > 1 and self.constant != 0:
            raise ValueError("constant must be zero when d + D > 1")
        backend = self.backend.lower()
        if backend not in {"auto", "statsmodels", "numpy"}:
            raise ValueError("backend must be one of 'auto', 'statsmodels', or 'numpy'")
        object.__setattr__(self, "backend", backend)

    @property
    def p(self) -> int:
        """Non-seasonal autoregressive order."""
        return self.order[0]

    @property
    def d(self) -> int:
        """Ordinary difference order."""
        return self.order[1]

    @property
    def q(self) -> int:
        """Non-seasonal moving-average order."""
        return self.order[2]

    @property
    def P(self) -> int:
        """Seasonal autoregressive order."""
        return self.seasonal_order[0]

    @property
    def D(self) -> int:
        """Seasonal difference order."""
        return self.seasonal_order[1]

    @property
    def Q(self) -> int:
        """Seasonal moving-average order."""
        return self.seasonal_order[2]

    def ar_polynomial(self, *, include_differences: bool = True) -> np.ndarray:
        """Return the AR-side lag polynomial in increasing lag order.

        Parameters
        ----------
        include_differences : bool, default True
            Include ordinary and seasonal difference polynomials. Set to false
            to obtain only the stationary AR factors used during simulation.

        Returns
        -------
        numpy.ndarray
            Polynomial coefficients beginning with one, with negative signs
            for the supplied AR terms. Seasonal cross-products are expanded.
        """
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
        """Return the expanded MA lag polynomial, beginning with one.

        Coefficients are in increasing lag order and use positive signs for
        the supplied MA terms. Ordinary and seasonal factors are multiplied.
        """
        nonseasonal = np.concatenate(([1.0], self.theta))
        seasonal = _seasonal_polynomial(self.Theta, self.frequency, sign=1.0)
        return np.polynomial.polynomial.polymul(nonseasonal, seasonal)

    def simulate(
        self,
        n: int = 100,
        *,
        n_start: int | None = None,
        rng: SeedLike = None,
        innov: Sequence[float] | np.ndarray | None = None,
        start_innov: Sequence[float] | np.ndarray | None = None,
    ) -> np.ndarray:
        """Simulate SARIMA, with scaled burn-in followed by integration.

        With n_start=None, an adaptive impulse-response energy check selects
        burn-in, including repeated and seasonal poles. Explicit zero starts
        from zero states. innov overrides n and supplies already-scaled
        output innovations; start_innov supplies burn-in innovations and its
        length determines n_start when that argument is omitted.

        Parameters
        ----------
        n : int, default 100
            Non-negative output length, excluding burn-in. Overridden by
            ``len(innov)`` when output innovations are supplied.
        n_start : int, optional
            Non-negative burn-in length. ``None`` selects an automatic length;
            zero starts from zero states without a burn-in period.
        rng : int or numpy.random.Generator, optional
            Seed or generator for any innovations not supplied explicitly.
            Existing generators advance in place; ``None`` creates a new one.
        innov : array_like, optional
            Finite output-period innovations, already multiplied by the desired
            noise scale. They are not multiplied by ``sigma`` again.
        start_innov : array_like, optional
            Finite, already-scaled burn-in innovations. Must match ``n_start``
            if both are supplied; otherwise their length determines it.

        Returns
        -------
        numpy.ndarray
            One-dimensional simulated path, including integration and the
            specified mean or drift. Zero output length returns an empty array.

        Raises
        ------
        ValueError
            If lengths or innovations are invalid, the AR part is nonstationary,
            or automatic burn-in fails within a window of ``2**22`` steps.
            Supply explicit ``n_start`` in the last case.
        ImportError
            If a required statsmodels component is unavailable.

        Notes
        -----
        The automatic search doubles impulse-response windows, stopping when
        the latter half contains at most ``1e-8`` of the window's energy. The
        returned burn-in is half that window. Its initial size covers the AR/MA
        orders and twelve time constants of the slowest pole, so short tails
        near AR-MA cancellation are not mistaken for convergence. This remains
        a numerical tail estimate, not an exact stationary initialization.

        Examples
        --------
        >>> model = arima_model(p=1, d=0, q=0, phi=[0.5], constant=0, sigma=1)
        >>> model.simulate(innov=[1, 0, 0], n_start=0).tolist()
        [1.0, 0.5, 0.25]
        """
        return self._simulate(n, n_start=n_start, rng=rng, innov=innov, start_innov=start_innov)

    def _prepare_simulator(self):
        if self.backend == "numpy":
            return self.simulate
        prepared = self._prepare_statsmodels()
        default_burnin = None

        def simulate_prepared(n=100, **kwargs):
            nonlocal default_burnin
            if (n and default_burnin is None and kwargs.get("n_start") is None
                    and kwargs.get("start_innov") is None):
                default_burnin = self._default_burnin()
            return self._simulate(n, prepared=prepared, default_burnin=default_burnin, **kwargs)

        return simulate_prepared

    def _simulate(self, n=100, *, n_start=None, rng=None, innov=None,
                  start_innov=None, prepared=None, default_burnin=None):
        n = _validate_order("n", n)
        if innov is not None:
            innov = self._validate_innovations(innov, "innov")
            n = innov.size
        if n == 0:
            return np.array([], dtype=float)
        _validate_sigma(self.sigma)
        for name, order in (("phi", self.p), ("theta", self.q), ("Phi", self.P), ("Theta", self.Q)):
            _validate_coefficients(name, getattr(self, name), order)
        if start_innov is not None:
            start_innov = self._validate_innovations(start_innov, "start_innov")
            if n_start is None:
                n_start = start_innov.size
        if n_start is None:
            n_start = self._default_burnin() if default_burnin is None else default_burnin
        else:
            n_start = _validate_order("n_start", n_start)
            self._check_stationarity()
        if start_innov is not None and start_innov.size != n_start:
            raise ValueError("start_innov must have length n_start")
        generator = get_rng(rng)
        if start_innov is None:
            start_innov = generator.normal(0.0, self.sigma, size=n_start)
        if innov is None:
            innov = generator.normal(0.0, self.sigma, size=n)
        innovations = np.concatenate((start_innov, innov))

        if self.backend in {"auto", "statsmodels"}:
            series = self._simulate_statsmodels(n=n, innovations=innovations, prepared=prepared)
        else:
            series = self._simulate_numpy(n=n, innovations=innovations)
        return self._integrate(series)

    @staticmethod
    def _validate_innovations(values: Sequence[float] | np.ndarray, name: str) -> np.ndarray:
        values = _as_1d(values)
        if np.any(~np.isfinite(values)):
            raise ValueError(f"{name} must be a finite numeric vector")
        return values

    def _check_stationarity(self) -> float:
        decay = np.inf
        for coefficients, period in ((self.phi, 1), (self.Phi, self.frequency)):
            if np.any(~np.isfinite(coefficients)):
                raise ValueError("AR coefficients must be finite")
            if not np.any(coefficients):
                continue
            radius = float(np.max(np.abs(np.roots(np.r_[1.0, -coefficients]))))
            if radius >= 1:
                raise ValueError("AR part of model must be stationary")
            if radius > 0:
                decay = min(decay, -np.log(radius) / period)
        return decay

    def _default_burnin(self) -> int:
        from statsmodels.tsa.arima_process import arma_impulse_response

        decay = self._check_stationarity()
        ar = np.trim_zeros(self.ar_polynomial(include_differences=False), "b")
        ma = np.trim_zeros(self.ma_polynomial(), "b")
        orders = ar.size + ma.size - 2
        if ar.size == 1:
            return int(orders)
        # Doubling windows cover delayed seasonal responses and repeated-pole
        # transients; normalize before squaring to avoid energy overflow.
        minimum = max(64, 2 * (orders + 1), int(np.ceil(12 / decay)))
        size = 1 << (minimum - 1).bit_length()
        while size <= 2**22:
            response = arma_impulse_response(ar, ma, leads=size)
            if not np.all(np.isfinite(response)):
                break
            energy = (response / np.max(np.abs(response))) ** 2
            if np.sum(energy[size // 2:]) <= 1e-8 * np.sum(energy):
                return size // 2
            size *= 2
        raise ValueError("automatic burn-in did not converge; supply n_start explicitly")

    def _integrate(self, series: np.ndarray) -> np.ndarray:
        for _ in range(self.d):
            series = np.cumsum(series)
        for _ in range(self.D):
            for offset in range(self.frequency):
                series[offset::self.frequency] = np.cumsum(series[offset::self.frequency])
        if self.d + self.D == 0:
            series = series + self.constant
        elif self.d + self.D == 1:
            series = series + self.constant * np.arange(1, series.size + 1)
        return series

    def _prepare_statsmodels(self):
        simulator, error = _load_statsmodels_arma()
        if simulator is None:
            raise ImportError("statsmodels is required for ARIMA simulation") from error
        return simulator, self.ar_polynomial(include_differences=False), self.ma_polynomial()

    def _simulate_statsmodels(self, *, n, innovations, prepared=None):
        """Filter supplied innovations using statsmodels, without dense states."""
        simulator, ar, ma = self._prepare_statsmodels() if prepared is None else prepared
        simulated = simulator(
            ar, ma, nsample=innovations.size,
            distrvs=lambda size: innovations,
        )
        return np.asarray(simulated, dtype=float).reshape(-1)[-n:]

    def _simulate_numpy(
        self,
        n: int,
        *,
        innovations: np.ndarray,
    ) -> np.ndarray:
        """Explicit NumPy simulator for the stationary ARMA part."""
        ar_poly = self.ar_polynomial(include_differences=False)
        ma_poly = self.ma_polynomial()
        ar_lag = ar_poly.size - 1
        ma_lag = ma_poly.size - 1
        start = max(ar_lag, ma_lag)
        innovations = np.r_[np.zeros(start), innovations]
        total = innovations.size
        y = np.zeros(total, dtype=float)

        for t in range(start, total):
            ar_part = 0.0
            if ar_lag:
                ar_part = -float(np.dot(ar_poly[1:], y[t - np.arange(1, ar_lag + 1)]))
            ma_part = 0.0
            if ma_lag:
                ma_part = float(np.dot(ma_poly[1:], innovations[t - np.arange(1, ma_lag + 1)]))
            y[t] = ar_part + innovations[t] + ma_part

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
    """Construct a Gaussian SARIMA model, sampling omitted parameters.

    Parameters
    ----------
    frequency : int, default 1
        Positive seasonal period, for example 12 for monthly data.
    p, q : int, optional
        Non-negative AR and MA orders. Inferred from supplied coefficient
        arrays; otherwise sampled from ``{0, 1, 2, 3}``.
    d : int, optional
        Ordinary difference order, sampled from ``{0, 1, 2}`` when omitted.
    P, Q : int, optional
        Seasonal AR and MA orders, inferred from supplied coefficients or
        sampled from ``{0, 1, 2}``. Zero for non-seasonal models.
    D : int, optional
        Seasonal difference order, sampled from ``{0, 1}`` subject to
        ``d + D <= 2``. Zero for non-seasonal models.
    constant : float, optional
        Mean with no differences or linear drift slope with one difference.
        Sampled uniformly from ``(-3, 3)`` unless there are two differences,
        when it must be zero. It is not an AR recursion intercept.
    phi, theta, Phi, Theta : array_like, optional
        Finite AR/MA coefficient vectors with lengths ``p``, ``q``, ``P`` and
        ``Q``. Missing AR and MA coefficients are sampled to be stationary and
        invertible, respectively. MA coefficients use positive polynomial signs.
    sigma : float, optional
        Innovation standard deviation, sampled uniformly from ``(1, 5)``.
        Its square must lie in the positive normal floating-point range.
    backend : {"auto", "statsmodels", "numpy"}, default "auto"
        ``auto`` and ``statsmodels`` use statsmodels' ARMA filter. ``numpy``
        uses an explicit recursion; there is no silent backend fallback.
    rng : int or numpy.random.Generator, optional
        Seed or generator used to sample missing parameters, not to store a
        future simulation seed. Existing generators advance in place; ``None``
        creates a new independently seeded generator.

    Returns
    -------
    ARIMAModel
        Model for ``simulate`` or ``generate``. Constructor and simulation RNGs
        are controlled separately unless the same generator is passed to both.

    Raises
    ------
    ValueError
        If orders, dimensions, differencing, noise scale or constants are
        invalid. Seasonal arguments require ``frequency > 1``.

    See Also
    --------
    ARIMAModel.simulate : Initialization and controlled-innovation semantics.
    stationary_ar : Sampling distribution for stationary AR coefficients.

    Examples
    --------
    >>> model = arima_model(p=1, d=0, q=0, phi=[0.5], constant=0, sigma=1)
    >>> model.simulate(12, rng=1).shape
    (12,)
    """
    generator = get_rng(rng)
    frequency = _validate_order("frequency", frequency)
    if frequency < 1:
        raise ValueError("frequency must be a positive integer")

    if sigma is None:
        sigma = float(generator.uniform(1.0, 5.0))
    _validate_sigma(sigma)

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
        if any(_validate_order("seasonal order", value) != 0 for value in seasonal_orders) or supplied_seasonal_coefficients:
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
    if not np.isfinite(constant):
        raise ValueError("constant must be finite")
    if d + D > 1 and constant != 0:
        raise ValueError("constant must be zero when d + D > 1")

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
