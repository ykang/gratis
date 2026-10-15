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
metadata. ETS uses only statsmodels; `backend="auto"` is an alias for that
backend and never falls back silently. ARIMA also defaults to statsmodels and
retains an explicit `backend="numpy"` simulator for debugging.

The statistical conventions match the R package:

- ARIMA `constant` is the mean with no differences, the slope of a linear
  drift with one difference, and zero with two differences. Burn-in applies
  to the stationary ARMA part before integration. R and Python use the same
  adaptive burn-in based on impulse-response energy, including repeated and
  seasonal poles. Doubling windows start at a power of two of at least 64 and
  twice the expanded AR plus MA orders plus one, and at least twelve time
  constants of the slowest pole to avoid false convergence near cancellation.
  Burn-in is half the first
  window whose latter half contains at most `1e-8` of its energy. This is an
  estimated tail criterion, not exact stationary initialization. Searches
  exceeding `2**22` steps fail and request explicit burn-in. All innovations
  use the requested `sigma`. Set `n_start=0` (R: `n.start=0`) for zero states;
  `innov` and `start_innov` (R: `start.innov`) allow controlled path comparisons.
  `generate()` calculates automatic burn-in once per batch in both languages.
- ARIMA and ETS `sigma` is an innovation standard deviation. Both languages
  reject scales whose squared variance is outside the positive normal
  floating-point range (approximately `1.49e-154 <= sigma <= 1.34e154`). Numeric MAR
  `sigmas` are innovation variances. R GARCH specifications supply standard
  deviations, which are squared before use by the MAR simulator. These R
  volatility paths are generated externally using latent GARCH shocks, not the
  resulting MAR residuals; this is not a coupled MAR-GARCH model.
- Stationary MAR components do not guarantee a stable mixture. Integrated
  components are allowed. MAR starts from zero states with a default 100-step
  burn-in; increase `n_start` (R: `n.start`) for persistent models, and assess
  mixture stability separately when a stationary distribution is required.
- MA coefficients use the positive-sign lag polynomial `1 + theta(B)`.
  Finite AR expansions preserve every seasonal lag. MA approximations warn
  if their tail has not decayed; increase `n_terms` when that occurs.
- ETS `season` runs from oldest to newest. Multiplicative seasonal factors
  must be positive; randomly generated factors are normalized to mean one.
- ETS uses the usual smoothing bounds and admissibility checks, with bounded
  resampling of omitted parameters. Multiplicative innovations are clipped
  below at -1, matching R's forecast simulator.
- Random stationary AR sampling above order two allows every possible count
  of real poles. It has full support, but is not uniform over coefficient space;
  seeded draws differ from the former fixed-root-count sampler.

Python MAR parameter arrays are copied and read-only. Construct a new model
instead of mutating coefficients; caller arrays cannot change an existing model.
Zero-order mixtures use vectorized simulation, and single-component MAR models
use a compiled AR filter. Both retain the original order of random draws.

Target search uses sorted breakpoints to cover the full mixture-weight simplex,
including zero-weight components. Both languages sample-standardize finite,
nonconstant paths using overflow-safe scaling; exact constants are unchanged.
An unmet target tolerance returns the best path with a warning reporting its
mean absolute feature error. Callback exceptions and incompatible feature
lengths raise errors; nonfinite candidate paths or features are rejected.
Mixture-normal covariances in both languages must be finite, symmetric and
positive definite. Sampling and density evaluation normalize nonnegative weights;
log densities are evaluated without first computing ordinary densities.

Python target search is serial by default. `simulate_target(..., parallel=True)`
evaluates candidates in threads. Callbacks must be thread-safe; speedups depend
on their releasing the GIL, and lightweight callbacks can be slower in threads.
Independent candidate random streams make serial and threaded results identical
for the same seed and a deterministic callback. This changes seeded target-search
results from versions that used one shared stream. A later invalid generation
does not discard an earlier valid best series. R also accepts caller-owned
clusters through `parallel`; callers must stop these clusters themselves.
R target search selects from paired final scores and series without retaining
unused generation histories or reevaluating feature callbacks.

`generate()` prepares statsmodels ARIMA/ETS backends once per batch, not once
per series. Preparation is not cached between calls, so array-valued model
parameters cannot leave a stale persistent backend. MAR recursion uses nonzero
lags, and scalar mixture simulation validates inputs once per path.
Conditional mixture densities use vectorized finite filters and skip inactive
components before calculating means. The public `ARIMAModel` and `ETSModel`
classes validate supplied coefficients, orders, scales, states, and applicable
smoothing constraints, just as the factory functions do.

R and NumPy have different random generators, so
the same seed does not produce identical paths across languages. Tests compare
fixed coefficients, mean/drift effects, and ETS paths with identical innovations.

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
package dependency, that means `statsmodels.tsa.arima_process.arma_generate_sample`.
It filters the expanded ARMA polynomials without allocating dense seasonal
state-space matrices. The same scaled innovations, burn-in, integration, and
mean/drift conventions are retained.
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

`ets_model()` uses `statsmodels.tsa.exponential_smoothing.ets.ETSModel`.
`backend="auto"` and `backend="statsmodels"` select the same implementation.
Missing or broken statsmodels installations raise an error:

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

Cross-language tests run when `Rscript` and the R packages `forecast` and
`polynom` are available. CI runs them alongside the Python tests.
The DGP regressions check repeated-pole initialization, supported root regimes,
and Monte Carlo MAR means, variances and lag-one correlations against analytic
values. These cover fixed reference models, not every possible generated DGP.

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
