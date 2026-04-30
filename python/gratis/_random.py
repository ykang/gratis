"""Random-number utilities shared by gratis modules."""

from __future__ import annotations

from typing import Optional, Union

import numpy as np

SeedLike = Optional[Union[int, np.random.Generator]]


def get_rng(seed: SeedLike = None) -> np.random.Generator:
    """Return a NumPy random generator from a seed or existing generator."""
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)
