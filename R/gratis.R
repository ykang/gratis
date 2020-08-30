#' Time Series Generation
#'
#' The tsgeneration package generates time series data based on MAR models.
#'
#' @docType package
#' @name tsgeneration
#' @importFrom stats as.ts coef dlnorm plnorm qlnorm rlnorm rmultinom rnorm runif stl ts
#' @importFrom stats diffinv fivenum frequency na.exclude na.omit sd window
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
#' @import shiny
NULL
# > NULL


utils::globalVariables(c(".","x","Season2"))
