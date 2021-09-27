#' Generating time series with controllable features.
#'
# Simulate one time series of length `length` with target features.
# feature_function is a function that returns a vector of features from a time series
# target is a vector of features of same length as that returned by feature_function
#' @param length length of the time series to be generated.
#' @param seasonal_periods Either a scalar or a numeric vector of length k containing
#' the number of seasonal periods for each component.
#' @param feature_function a function that returns a vector of features from a time series
#' @param target target feature values of the same length as that returned by \code{feature_function()}
#' @param parallel An optional argument which allows to specify if the Genetic Algorithm
#'     should be run sequentially or in parallel.
#' @return A time-series object of class "ts" or "msts".
#' @author Yanfei Kang and Rob J Hyndman
#' @examples
#' \dontrun{
#' my_features <- function(y) {
#'   c(tsfeatures::entropy(y), acf = acf(y, plot = FALSE)$acf[2:3, 1, 1])
#' }
#' y <- simulate_target(
#'   length = 60, feature_function = my_features,
#'   target = c(0.5, 0.9, 0.8)
#' )
#' my_features(y)
#' plot(y)}
#' @export
#'
simulate_target <- function(length, seasonal_periods = 1, feature_function, target, parallel = FALSE) {
  # Set min and max
  #
  if (length(seasonal_periods == 1)) {
    ga_length <- 24
  } else {
    ga_length <- 48
  }
  ga_min <- rep(-3, ga_length)
  ga_max <- rep(3, ga_length)
  # ga_min = 0
  # ga_max = 1

  # Test lengths
  feature_length <- length(feature_function(msts(rnorm(100), seasonal.periods = seasonal_periods)))
  if (length(target) != feature_length) {
    stop(paste("target should be of length", feature_length))
  }

  GA <- ga_ts(
    type = "real-valued", fitness = fitness_mar,
    length = length, seasonal_periods, ncomponents = 3,
    feature_function = feature_function, target = target,
    n = length,
    min = ga_min,
    max = ga_max,
    parallel = parallel, popSize = 50, maxiter = 200,
    pmutation = 0.3, pcrossover = 0.8, maxFitness = 0,
    run = 30, keepBest = TRUE, monitor = NULL
  )
  forecast::msts(as.matrix(tail(GA@bestSol, 1)[[1]])[, 1], seasonal.periods = seasonal_periods)
}

fitness_mar <- function(pars, length, seasonal_periods, ncomponents, feature_function, target) {
  npar <- length(pars)
  k <- ncomponents
  p <- npar / k - 4
  ar <- matrix(pars[seq((p + 1) * k)], ncol = k)
  sigmas <- abs(pars[(p + 1) * k + seq(k)])
  weights <- abs(pars[(p + 2) * k + seq(k)])
  x <- simulate(mar_model(ar = ar, sigmas = sigmas, weights = weights, seasonal_periods = seasonal_periods),
    nsim = length
  )
  return(list(
    value = -sum((feature_function(x) - target)^2),
    x = x
  ))
}
