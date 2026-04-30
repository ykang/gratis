"""Common generation entry points."""

from __future__ import annotations

import numpy as np

from ._random import SeedLike, get_rng


def simulate(model, n: int = 100, *, nsim: int | None = None, **kwargs) -> np.ndarray:
    """Simulate one series from a model object with a `simulate` method."""
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
    """Generate multiple synthetic series as a length x nseries array."""
    if length < 0:
        raise ValueError("length must be non-negative")
    if nseries < 1:
        raise ValueError("nseries must be positive")

    generator = get_rng(rng)
    out = np.empty((length, nseries), dtype=float)
    for index in range(nseries):
        out[:, index] = simulate(model, n=length, rng=generator, **kwargs)
    return out
