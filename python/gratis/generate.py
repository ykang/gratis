"""Common generation entry points."""

from __future__ import annotations

from functools import partial

import numpy as np

from ._random import SeedLike, get_rng
from .arima import _validate_order


def simulate(model, n: int = 100, *, nsim: int | None = None, **kwargs) -> np.ndarray:
    """Simulate one path from a specified model.

    Parameters
    ----------
    model : MARModel, ARIMAModel, ETSModel or compatible object
        Object providing a ``simulate(n=..., **kwargs)`` method.
    n : int, default 100
        Number of observations, excluding any model-specific burn-in.
    nsim : int, optional
        R-compatible alias for ``n``. When supplied, it takes precedence.
    **kwargs
        Forwarded to the model's simulation method, including ``rng`` and any
        supported initialization controls. Not every model accepts the same
        arguments; for example, ETS has no burn-in argument.

    Returns
    -------
    numpy.ndarray
        One-dimensional floating-point series for the built-in models. No R
        time-series or calendar metadata is attached.

    Raises
    ------
    NotImplementedError
        If the object does not provide a simulation method.

    See Also
    --------
    generate : Generate several paths as columns of an array.

    Examples
    --------
    >>> from gratis import arima_model, simulate
    >>> model = arima_model(p=1, d=0, q=0, phi=[0.5], constant=0, sigma=1)
    >>> simulate(model, n=12, rng=1).shape
    (12,)
    """
    if nsim is not None:
        n = nsim
    method = getattr(model, "simulate", None)
    if method is None:
        raise NotImplementedError(f"{type(model).__name__} does not provide simulate()")
    return np.asarray(method(n=n, **kwargs), dtype=float)


def generate(
    model,
    length: int = 100,
    nseries: int = 10,
    *,
    rng: SeedLike = None,
    **kwargs,
) -> np.ndarray:
    """Generate a batch of paths, with one series per column.

    Parameters
    ----------
    model : MARModel, ARIMAModel, ETSModel or compatible object
        Model from which each path is simulated.
    length : int, default 100
        Non-negative number of observations per path, excluding burn-in.
    nseries : int, default 10
        Positive number of paths to generate.
    rng : int or numpy.random.Generator, optional
        Integer seed or existing generator. A seed initializes one generator
        for the whole batch; successive paths consume successive draws. An
        existing generator is used and advanced in place. ``None`` creates a
        new generator and does not use NumPy's legacy global RNG.
    **kwargs
        Arguments forwarded to each simulation, such as MAR ``n_start`` or
        ARIMA ``innov`` and ``start_innov``. Supplied innovations must produce
        exactly ``length`` observations.

    Returns
    -------
    numpy.ndarray
        Floating-point array of shape ``(length, nseries)``, including when
        ``nseries=1``. Paths do not carry calendar or seasonal metadata.

    Raises
    ------
    ValueError
        If the dimensions are invalid or a simulation has an incompatible
        length. Model-specific errors propagate.

    Notes
    -----
    Statsmodels preparation is reused within a batch, not across batches.
    Automatic ARIMA burn-in is also reused for the statsmodels backend.

    Examples
    --------
    >>> from gratis import arima_model, generate
    >>> model = arima_model(p=1, d=0, q=0, phi=[0.5], constant=0, sigma=1)
    >>> generate(model, length=12, nseries=3, rng=1).shape
    (12, 3)
    """
    length = _validate_order("length", length)
    nseries = _validate_order("nseries", nseries)
    if nseries < 1:
        raise ValueError("nseries must be positive")

    generator = get_rng(rng)
    out = np.empty((length, nseries), dtype=float)
    prepare = getattr(model, "_prepare_simulator", None)
    method = prepare() if prepare is not None and length else partial(simulate, model)
    for index in range(nseries):
        out[:, index] = method(n=length, rng=generator, **kwargs)
    return out
