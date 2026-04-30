#' Generate time series with target features
#'
#' This function is deprecated. Use \code{\link{generate_target}()} instead.
#'
#' @param n Number of time series to generate.
#' @param ts.length Length of each generated time series.
#' @param freq Seasonal frequency, or seasonal periods for
#'   multiple-seasonal data.
#' @param seasonal Integer seasonality flag: 0 for non-seasonal data, 1 for
#'   single-seasonal data, and 2 for multiple-seasonal data.
#' @param features Character vector of feature function names passed to
#'   \code{\link[tsfeatures]{tsfeatures}()}.
#' @param selected.features Character vector naming the features to match.
#' @param target Numeric vector of target feature values in the same order as
#'   \code{selected.features}.
#' @param parallel Logical value or parallel backend specification indicating
#'   whether the genetic algorithm should evaluate candidates in parallel.
#' @param output_format Output format: \code{"list"} for a \code{ts} or
#'   \code{msts} object, or \code{"tsibble"} for a tsibble.
#' @return Generated time series in the requested \code{output_format}.
#' @author Yanfei Kang
#' @examples
#' \dontrun{
#' library(tsfeatures)
#' x <- generate_ts_with_target(
#'   n = 1, ts.length = 60, freq = 1, seasonal = 0,
#'   features = c("entropy", "stl_features"), selected.features = c("entropy", "trend"),
#'   target = c(0.6, 0.9), parallel = FALSE
#' )
#' forecast::autoplot(x)
#' }
#' @export
generate_ts_with_target <- function(n, ts.length, freq, seasonal, features, selected.features, target, parallel = TRUE, output_format = "list") {
  if (!requireNamespace("tsfeatures", quietly = TRUE)) {
    stop("Package 'tsfeatures' is required to generate time series with target features.", call. = FALSE)
  }

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
      type = "real-valued", fitness = fitness_ts, features = features, seasonal = seasonal,
      ts.length, freq, target, 3, selected.features,
      n = ts.length,
      min = ga_min,
      max = ga_max,
      parallel = parallel, popSize = 30, maxiter = 100,
      pmutation = 0.3, pcrossover = 0.8, maxFitness = -sqrt(0.01 * length(selected.features)),
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
  # New content
  output <- if (output_format == "list") {
    evolved.ts
  } else if (output_format == "tsibble") {
    tsibble::as_tsibble(evolved.ts)
  }
  return(output)
}
