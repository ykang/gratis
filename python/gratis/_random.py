"""Random-number utilities shared by gratis modules."""

from __future__ import annotations

from typing import Optional, Union

import numpy as np

SeedLike = Optional[Union[int, np.random.Generator]]


def get_rng(seed: SeedLike = None) -> np.random.Generator:
    """Return an existing generator or create one with ``default_rng``.

    An existing generator is returned unchanged, so callers advance its state.
    Integer seeds reproduce draws within the Python implementation. ``None``
    creates an independently seeded generator; NumPy's legacy global RNG is
    neither consulted nor modified. Equal seeds do not imply equal R draws.
    """
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)
