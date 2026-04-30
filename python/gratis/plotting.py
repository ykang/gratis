"""Plotting helpers for generated time series."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from .generate import generate


def _series_matrix(data: Sequence[float] | np.ndarray) -> np.ndarray:
    values = np.asarray(data, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    elif values.ndim != 2:
        raise ValueError("data must be a one-dimensional series or a two-dimensional matrix")
    if values.shape[0] == 0:
        raise ValueError("data must contain at least one observation")
    return values


def _plot_labels(labels: Sequence[str] | None, nseries: int) -> list[str]:
    if labels is None:
        return [f"Series {index + 1}" for index in range(nseries)]
    labels = list(labels)
    if len(labels) != nseries:
        raise ValueError("labels must have one entry for each plotted series")
    return labels


def _plot_index(index: Sequence[Any] | np.ndarray | None, length: int) -> np.ndarray:
    if index is None:
        return np.arange(length)
    index_values = np.asarray(index)
    if index_values.shape[0] != length:
        raise ValueError("index length must match the number of observations")
    return index_values


def _get_axes(ax: Any | None):
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "plotting requires matplotlib. Install it with "
            "`python -m pip install -e './python[plot]'`."
        ) from exc
    _, ax = plt.subplots()
    return ax


def plot_series(
    data: Sequence[float] | np.ndarray,
    *,
    ax: Any | None = None,
    index: Sequence[Any] | np.ndarray | None = None,
    labels: Sequence[str] | None = None,
    title: str | None = None,
    max_series: int | None = None,
    legend: bool = True,
    alpha: float = 0.9,
    linewidth: float = 1.5,
    **plot_kwargs: Any,
):
    """Plot one or more generated time series.

    Parameters
    ----------
    data
        A one-dimensional series or a two-dimensional array shaped
        ``length x nseries``.
    ax
        Optional matplotlib-like axes. When omitted, matplotlib is imported
        lazily and a new axes is created.
    index
        Optional x-axis values with one entry per observation.
    labels
        Optional labels with one entry per plotted series.
    title
        Optional axes title.
    max_series
        Plot at most this many columns from a multi-series matrix.
    legend
        Whether to draw a legend for multi-series data.
    alpha, linewidth, **plot_kwargs
        Passed to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes used for plotting.
    """
    values = _series_matrix(data)
    if max_series is not None:
        max_series = int(max_series)
        if max_series < 1:
            raise ValueError("max_series must be positive")
        values = values[:, :max_series]

    x = _plot_index(index, values.shape[0])
    series_labels = _plot_labels(labels, values.shape[1])
    ax = _get_axes(ax)

    for column, label in enumerate(series_labels):
        ax.plot(
            x,
            values[:, column],
            label=label,
            alpha=alpha,
            linewidth=linewidth,
            **plot_kwargs,
        )

    if title is not None:
        ax.set_title(title)
    ax.set_xlabel("Index")
    ax.set_ylabel("Value")
    if legend and values.shape[1] > 1:
        ax.legend()
    return ax


def plot_generated(
    model,
    *,
    length: int = 100,
    nseries: int = 1,
    ax: Any | None = None,
    rng=None,
    **kwargs: Any,
):
    """Generate series from a model and plot them."""
    plot_options = {
        key: kwargs.pop(key)
        for key in list(kwargs)
        if key in {"index", "labels", "title", "max_series", "legend", "alpha", "linewidth"}
    }
    data = generate(model, length=length, nseries=nseries, rng=rng, **kwargs)
    return plot_series(data, ax=ax, **plot_options)
