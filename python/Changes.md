# gratis Python 1.9

## Major Bug Fixes

* Aligned ARIMA mean/drift conventions, MA coefficient signs, and seasonal
  polynomial expansion with R. Preserved long seasonal lags and zero-sum
  seasonal coefficients; unconverged MA approximations now warn.
* Corrected ARIMA innovation scaling during startup. Automatic burn-in uses an
  adaptive impulse-energy estimate for repeated poles and near AR-MA
  cancellation. Nonconvergence raises an error requesting explicit burn-in;
  explicit burn-in bypasses the search without bypassing stationarity checks.
* Random stationary AR models above order two now allow every possible count
  of real roots. Sampling has full support but is not uniform over coefficient
  space.
* Added bounded ETS parameter resampling and admissibility checks, aligned
  seasonal-state ordering with R, and required positive multiplicative states
  and seasonal factors.
* Target-search weights cover the full mixture simplex, including zero weights.
  Unmet tolerances return the best series with a warning; callback errors
  propagate and nonfinite candidates are rejected. A valid incumbent is retained
  if a later generation is entirely invalid.
* Made sample standardization overflow-safe, preserving exactly constant series
  and small variations around large offsets, consistently with R.
* Mixture covariance matrices must be finite, symmetric and positive definite.
  Weights are normalized and log densities remain stable in distribution tails.
* Public ARIMA and ETS classes validate supplied parameters directly, including
  coefficients, orders, noise scales, states and applicable smoothing bounds.

## Performance

* Added sparse MAR lag evaluation, vectorized zero-order mixtures, and compiled
  single-component AR filtering while retaining the order of random draws.
* Vectorized conditional mixture densities, skipped zero-weight conditional
  means, and avoided repeated scalar-mixture validation and variance allocation.
* ARIMA uses statsmodels' ARMA filter without dense state-space setup.
  Statsmodels ARIMA/ETS preparation and automatic ARIMA burn-in are reused within
  generation batches, rather than recalculated for each path.
* Target search supports threaded candidate evaluation with reproducible random
  streams. Serial and threaded results match for the same seed and a deterministic,
  thread-safe feature callback; speedups depend on the callback workload.

## Compatibility

* MAR parameter arrays are copied and read-only; construct a new model to change
  parameters instead of mutating arrays.
* ARIMA/ETS noise scales whose squared variance lies outside the positive normal
  floating-point range are rejected, matching R's validation.
* Seeded random-model draws, default ARIMA paths, and target-search results can
  differ from earlier versions because sampling and initialization have changed.
* ETS uses statsmodels only. ARIMA also defaults to statsmodels and retains an
  explicit `backend="numpy"` alternative; neither silently falls back.
* Python sources and documentation remain excluded from the R package archive.

## Documentation And Tests

* Expanded NumPy-style docstrings with parameter types, array shapes, noise
  units, RNG behavior, initialization, warnings, cross-references and executable
  examples. Updated the Python README with statistical conventions and limits.
* Added pytest regressions for numerical safety, simulation boundaries, seeded
  recursion, parameter validation and ownership, and threaded target searches.
* Cross-language checks compare controlled paths and conventions. Monte Carlo
  checks compare reference MAR means, variances and lag-one correlations with
  analytic values.
