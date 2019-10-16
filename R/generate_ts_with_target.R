#' Generating time series with controllable features.
#'
#' @param n number of time series to be generated.
#' @param ts.length length of the time series to be generated.
#' @param freq frequency of the time series to be generated.
#' @param seasonal 0 for non-seasonal data, 1 for single-seasonal data, and 2 for multiple seasonal data.
#' @param features a vector of function names.
#' @param selected.features selected features to be controlled.
#' @param target target feature values.
#' @param parallel An optional argument which allows to specify if the Genetic Algorithm
#'     should be run sequentially or in parallel.
#' @return A time-series object of class "ts" or "msts".
#' @author Yanfei Kang
#' @examples
#' library(tsfeatures)
#' x <- generate_ts_with_target(n = 1, ts.length = 60, freq = 1, seasonal = 0,
#'   features = c('entropy', 'stl_features'), selected.features = c('entropy', 'trend'),
#'   target=c(0.6, 0.9),  parallel=FALSE)
#' forecast::autoplot(x)
#' @export
generate_ts_with_target <- function(n, ts.length, freq, seasonal, features, selected.features, target, parallel=TRUE) {
  ga_min <-
    if (seasonal == 0) {
      c(rep(0, 10))
    } else if (seasonal == 1) {
      c(rep(0, 17))
    } else {
      c(rep(0, 35))
    }
  ga_max <-
    if (seasonal == 0) {
      c(rep(1, 10))
    } else if (seasonal == 1) {
      c(rep(1, 17))
    } else {
      c(rep(1, 35))
    }
  evolved.ts <- c()
  while (ifelse(is.null(dim(evolved.ts)), 0 < 1, dim(evolved.ts)[2] < n)) {
    GA <- ga_ts(
      type = "real-valued", fitness = gratis::fitness_ts, features = features, seasonal = seasonal,
      ts.length, freq, target, 3, selected.features,
      n = ts.length,
      min = ga_min,
      max = ga_max,
      parallel = parallel, popSize = 30, maxiter = 100,
      pmutation = 0.3, pcrossover = 0.8, maxFitness = -sqrt(0.01*length(selected.features)),
      run = 30, keepBest = TRUE, monitor = GA::gaMonitor
    )
    evolved.ts.new <-
      unique(do.call(
        cbind,
        eval(parse(text = paste("list(", paste("GA@bestSol[[GA@iter - ", 0:(GA@run - 1), "]]", sep = "", collapse = ","), ")")))
      ), MARGIN = 2)
    evolved.ts <- cbind(evolved.ts, evolved.ts.new)
  }
  if (length(freq) == 1) {
    evolved.ts <- ts(evolved.ts[, 1:n], frequency = freq)
  } else {
    evolved.ts <- msts(evolved.ts[, 1:n], seasonal.periods = freq)
  }
  return(evolved.ts)
}
