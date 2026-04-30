# gratis Python

This directory contains the Python package for the `gratis` project. The R
package remains at the repository root, while Python packaging starts in this
subdirectory and source code lives in `python/gratis`.

The Python package is intentionally excluded from `R CMD build` by the root
`.Rbuildignore`, so the repository can be used as both an R package and a
Python package without shipping Python sources in the R package tarball.

## Status

The Python implementation currently provides the main simulation workflows:

- `gratis.mixture`: mixture-normal random generation and densities.
- `gratis.arima`: SARIMA coefficient utilities, ARIMA model construction, and
  ARIMA simulation backed by `statsmodels`.
- `gratis.mar`: mixture autoregressive model construction and simulation.
- `gratis.ets`: ETS model construction and simulation backed by `statsmodels`.
- `gratis.generate`: common `simulate()` and `generate()` helpers.
- `gratis.plotting`: plotting helpers for generated arrays and models.
- `gratis.target`: target-feature generation using an evolutionary search over
  MAR parameters.

Python outputs are NumPy arrays. They do not carry R `ts`, `msts`, or `tsibble`
metadata. ARIMA and ETS have local NumPy fallback simulators available through
`backend="numpy"` when useful for debugging.

## Installation

From the repository root:

```sh
python -m pip install -e ./python
```

For development and tests:

```sh
python -m pip install -e "./python[dev]"
```

You can also work directly from the repository without installing:

```sh
PYTHONPATH=python python -c "import gratis; print(gratis.__version__)"
```

## Quick Start

```python
import numpy as np
import gratis

model = gratis.mar_model(
    phi=np.array([[0.8, 0.6], [0.0, 0.3]]),
    d=0,
    sigmas=[1.0, 2.0],
    weights=[0.8, 0.2],
)

series = gratis.simulate(model, n=100, rng=1)
many = gratis.generate(model, length=100, nseries=5, rng=1)
```

Plot generated arrays or generate and plot in one step:

```python
gratis.plot_series(many, title="Generated MAR series")
gratis.plot_generated(model, length=100, nseries=5, rng=1)
```

## Model Examples

ARIMA:

```python
model = gratis.arima_model(
    p=1,
    d=0,
    q=1,
    phi=[0.4],
    theta=[0.2],
    constant=0.0,
    sigma=1.0,
)

series = model.simulate(n=100, rng=1)
```

`arima_model()` uses `backend="auto"` by default. Because `statsmodels` is a
package dependency, that normally means `statsmodels.tsa.statespace.sarimax.SARIMAX`.
To require the maintained backend and fail if it is unavailable:

```python
model = gratis.arima_model(
    p=1,
    d=0,
    q=1,
    phi=[0.4],
    theta=[0.2],
    constant=0.0,
    sigma=1.0,
    backend="statsmodels",
)
```

ETS:

```python
model = gratis.ets_model(
    error="A",
    trend="N",
    seasonal="N",
    alpha=0.2,
    level=1.0,
    sigma=1.0,
)

series = gratis.simulate(model, n=100, rng=1)
```

`ets_model()` uses `backend="auto"` by default. Because `statsmodels` is a
package dependency, that normally means
`statsmodels.tsa.exponential_smoothing.ets.ETSModel`. To require the maintained
backend and fail if it is unavailable:

```python
model = gratis.ets_model(
    error="A",
    trend="N",
    seasonal="N",
    alpha=0.2,
    level=1.0,
    sigma=1.0,
    backend="statsmodels",
)
```

Target-feature generation:

```python
def features(x):
    return np.array([np.mean(x), np.std(x)])

series = gratis.simulate_target(
    length=100,
    feature_function=features,
    target=np.array([0.0, 1.0]),
    rng=1,
)
```

## Tests

Run the Python tests from `python/`:

```sh
python -m pytest -q
```

Or from the repository root:

```sh
PYTHONPATH=python python -m pytest -q python/tests
```

The root R package should still exclude this directory:

```sh
R CMD build --no-build-vignettes .
tar -tzf gratis_*.tar.gz | rg '^gratis/python|python/'
```

The final `tar` command should print no matches.
