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

* Added new parameter 'output_format', give users an option to transform time series output data into a tsibble format.
* Update vignette content

## Bug fixes

* Corrected dependencies in R code, which is namespaces in Imports field not imported from: 'rlang' 'shinydashboard'

# gratis 0.2.0

* First release.

## New features

* The R package gratis (previously known as tsgeneration) provides efficient algorithms for generating time series with diverse and controllable characteristics.
