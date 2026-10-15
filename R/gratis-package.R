#' gratis: Generate synthetic time series
#'
#' Generate synthetic univariate time series with diverse or controllable
#' characteristics. The package supports mixture autoregressive (MAR), ARIMA,
#' and ETS data-generating processes, and includes a genetic-algorithm workflow
#' for searching MAR parameters that produce user-specified time-series
#' features.
#'
#' The main workflow is to construct a model with \code{\link{mar_model}()},
#' \code{\link{arima_model}()}, or \code{\link{ets_model}()}, then simulate one
#' path with \code{\link[stats]{simulate}()} or many paths with
#' \code{\link[generics]{generate}()}. To generate series with target feature
#' values, use \code{\link{simulate_target}()} or
#' \code{\link{generate_target}()}.
#'
#' ARIMA and ETS \code{sigma} parameters are innovation standard deviations;
#' numeric MAR \code{sigmas} are innovation variances. Specifying a stationary
#' component AR polynomial alone does not ensure a stable MAR mixture or an
#' adequate burn-in period. See each model's help for initialization controls.
#'
#' The Python implementation in this repository follows the same statistical
#' conventions, but returns NumPy arrays instead of R time-series objects.
#' Equal numeric seeds do not imply identical random draws across languages.
#'
#' @aliases gratis-package
#' @importFrom stats as.ts coef dlnorm plnorm qlnorm rlnorm rmultinom rnorm runif stl ts
#' @importFrom stats diffinv fivenum frequency na.exclude na.omit sd window simulate
#' @importFrom generics generate
#' @importFrom utils head write.csv str flush.console globalVariables
#' @importFrom graphics points polygon title
#' @importFrom grDevices adjustcolor
#' @importFrom methods new slotNames
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom polynom polynomial
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom forecast trendcycle msts autoplot
#' @importFrom dplyr mutate select bind_cols
#' @importFrom magrittr "%>%"
#' @importFrom foreach foreach "%dopar%" "%do%"
#' @importFrom GA ga gaControl gaMonitor
#' @importFrom GA startParallel garun
#' @importFrom doRNG "%dorng%"
#' @keywords internal

"_PACKAGE"

NULL
