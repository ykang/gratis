"""Generate synthetic time series from MAR, ARIMA and ETS models.

Construct a model with ``mar_model``, ``arima_model`` or ``ets_model``, then use
``simulate`` for one path or ``generate`` for a matrix of paths. ``plot_series``
and ``plot_generated`` return Matplotlib axes. Target-feature searches use
``simulate_target`` and ``generate_target``.

Outputs are NumPy arrays: observations occupy rows and separate series occupy
columns. Randomized entry points accept an integer seed or NumPy generator via
``rng``. R-compatible statistical conventions do not imply identical seeded
draws or identical optimization results across languages.
"""

from ._version import __version__
from .arima import ARIMAModel, arima_model, pi_coefficients, stationary_ar
from .ets import ETSModel, ets_admissible, ets_model
from .generate import generate, simulate
from .mar import MARModel, mar_model
from .mixture import dmixnorm, dmixnorm_ts, rmixnorm, rmixnorm_ts
from .plotting import plot_generated, plot_series
from .target import generate_target, scalets, simulate_target

__all__ = [
    "__version__",
    "ARIMAModel",
    "ETSModel",
    "MARModel",
    "arima_model",
    "dmixnorm",
    "dmixnorm_ts",
    "ets_admissible",
    "ets_model",
    "generate",
    "generate_target",
    "mar_model",
    "pi_coefficients",
    "plot_generated",
    "plot_series",
    "rmixnorm",
    "rmixnorm_ts",
    "scalets",
    "simulate",
    "simulate_target",
    "stationary_ar",
]
