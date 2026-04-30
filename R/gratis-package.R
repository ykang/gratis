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
#' @importFrom tsfeatures tsfeatures
#' @importFrom doRNG "%dorng%"
#' @importFrom fGarch garchSpec garchSim
#' @import shiny
#' @import feasts
#' @keywords internal

"_PACKAGE"

NULL
