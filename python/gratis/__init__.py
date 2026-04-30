"""Python implementation of the gratis time-series generation package."""

from ._version import __version__
from .arima import ARIMAModel, arima_model, pi_coefficients, stationary_ar
from .ets import ETSModel, ets_admissible, ets_model
from .generate import generate, simulate
from .mar import MARModel, mar_model
from .mixture import dmixnorm, dmixnorm_ts, rmixnorm, rmixnorm_ts
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
    "rmixnorm",
    "rmixnorm_ts",
    "scalets",
    "simulate",
    "simulate_target",
    "stationary_ar",
]
