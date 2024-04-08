#' gratis: GeneRAting TIme Series with diverse and controllable characteristics
#'
#' The gratis package generates synthetic time series data based on various
#' univariate time series models including MAR, ARIMA and ETS processes.
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
