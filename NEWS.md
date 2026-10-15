# gratis 1.9

## Major Bug Fixes

* Aligned R and Python ARIMA mean/drift conventions, MA coefficient signs, and
  seasonal polynomial expansion. Preserved long seasonal lags and zero-sum
  seasonal coefficients; unconverged MA approximations now warn.
* Corrected ARIMA innovation scaling during startup. Automatic burn-in now uses
  an adaptive impulse-energy estimate for repeated poles and near AR-MA
  cancellation. Nonconvergence raises an error requesting explicit burn-in;
  explicit burn-in bypasses the search without bypassing stationarity checks.
* Random stationary AR models above order two now allow every possible count
  of real roots instead of excluding valid stationary regimes. Sampling has
  full support but is not uniform over coefficient space.
* Fixed ETS parameter resampling that could hang on rejected parameters.
  Added bounded admissibility checks, corrected seasonal-state ordering, and
  required positive multiplicative states and seasonal factors.
* Target-search weights now cover the full mixture simplex, including zero
  weights. Unmet tolerances return the best series with a warning; callback
  errors propagate and nonfinite candidates are rejected. Python preserves a
  valid incumbent if a later generation is entirely invalid.
* Made sample standardization overflow-safe and consistent across languages,
  preserving exactly constant series and small variations around large offsets.
  Mixture covariance matrices must be finite, symmetric and positive definite;
  weights are normalized and log densities remain stable in distribution tails.
* Fixed R single-component target indexing, zero-burn-in MAR lengths, and short
  mixture-series indexing. Corrected variance scaling for external GARCH paths.
* Repaired R GA selection and crossover bookkeeping, scalar local optimization,
  parallel and optimized population updates, oversized suggestions, elitism,
  nonuniform mutation, and permutation ranges. Parameters, scores and stored
  series remain paired; supplied clusters remain caller-owned.
* R batch generation seeds once per batch, restores the previous RNG state,
  and preserves calendar indices in target-generated data.

## Performance

* Added sparse MAR lag evaluation, vectorized zero-order mixtures, and compiled
  single-component AR filtering while retaining the order of random draws.
* Vectorized conditional mixture densities, skipped zero-weight conditional
  means, and avoided repeated scalar-mixture validation and variance allocation.
* Python ARIMA uses statsmodels' ARMA filter without dense state-space setup.
  Python ARIMA/ETS backend preparation and R/Python automatic ARIMA burn-in are
  reused within each generation batch.
* R target search selects paired final scores and series without retaining
  unused generation histories or reevaluating final feature callbacks. Python
  supports threaded candidate evaluation with reproducible random streams.

## Compatibility

* Removed the bundled Shiny applications, `app_gratis()` launcher, and app-only
  suggested dependencies. Python remains excluded from the R package archive.
* Python MAR parameter arrays are copied and read-only; construct a new model
  to change parameters. Public Python ARIMA and ETS classes validate inputs.
* R and Python reject ARIMA/ETS noise scales whose squared variance lies outside
  the positive normal floating-point range, rather than silently overflowing
  or underflowing.
* Seeded random-model draws, default ARIMA paths, and Python target-search
  results can differ from earlier versions because sampling and initialization
  have changed. Python ETS uses statsmodels only; ARIMA retains an explicit
  `backend="numpy"` alternative.

## Documentation And Tests

* Updated roxygen documentation and both READMEs with parameter conventions,
  initialization limits, and the distinction between external GARCH volatility
  and residual-coupled MAR-GARCH dynamics. Updated `simulate_target()` examples
  to use `dplyr::reframe()` instead of `summarise()` for dplyr 1.2.0.
* Added testthat and pytest regressions for numerical safety, simulation
  boundaries, seeded recursion, GA bookkeeping, and real parallel workers.
  Cross-language checks compare controlled paths and conventions; Monte Carlo
  checks compare reference MAR means, variances and lag-one correlations with
  analytic values.

# gratis 1.0.7

* Fixed dependency issue on CRAN.

# gratis 1.0.5

* Fixed the 'docType' issue.

# gratis 1.0.3

* Fixed CRAN checking issues.

# gratis 1.0.0

* Improved documentation
* Added new functions for MAR simulation: mar_model(), simulate.mar()
* Added ets_model() and arima_model() to randomly specify ETS and ARIMA models
* Added generate.ets, generate.mar, and generate.Arima to return tsibbles
* Added generate_target() and simulate_target() functions.
* Deprecated generate_msts(), generate_ts(), and generate_ts_with_target().
* Removed some redundant functions

# gratis 0.2.1

## Improvements

* Added new parameter `output_format`, giving users an option to transform time series output data into a `tsibble` format.
* Update vignette content

## Bug fixes

* Corrected dependencies in R code, which is namespaces in Imports field not imported from: 'rlang' 'shinydashboard'

# gratis 0.2.0

* First release.

## New features

* The R package `gratis` (previously known as tsgeneration) provides efficient algorithms for generating time series with diverse and controllable characteristics.
