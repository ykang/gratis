#' Generating time series with controllable features using MAR models
#'
#' Simulate one time series of length `length` from a MAR model with target features.
#' \code{feature_function} is a function that returns a vector of features from a time series,
#' while \code{target} is a vector of features of the same length as that returned by
#' \code{feature_function}.
#' @param length length of the time series to be generated.
#' @param seasonal_periods Either a scalar or a numeric vector of length k containing
#' the number of seasonal periods for each component.
#' @param feature_function a function that returns a vector of features from a time series.
#' @param target target feature values of the same length as that returned by \code{feature_function()}.
#' @param model character string indicating what model should be used to generate the data.
#' @param p integer specifying order of autoregression to use for MAR models. Default is 2.
#' @param k integer specifying number of components to use for MAR models. Default is 3
#' unless there are multiple seasonal periods specified.
#' @param tolerance average tolerance per feature. The genetic algorithm will attempt
#' to find a solution where the average difference between the series features and target
#' is less than \code{tolerance}. A larger value will give a faster but less precise
#' solution.
#' @param parallel An optional argument which allows to specify if the Genetic Algorithm
#'     should be run sequentially or in parallel.
#' @return A time-series object of class "ts" or "msts".
#' @author Yanfei Kang and Rob J Hyndman
#' @examples
#' \dontrun{
#' set.seed(1)
#' library(tsfeatures)
#' my_features <- function(y) {
#'   c(entropy(y), acf = acf(y, plot = FALSE)$acf[2:3, 1, 1])
#' }
#' # Specify target features
#' y <- simulate_target(
#'   length = 60, p = 1, 
#'   feature_function = my_features, target = c(0.5, 0.9, 0.8)
#' )
#' my_features(y)
#' plot(y)
#' # Simulate a time series similar to an existing series
#' my_features <- function(y) {
#'   c(min = min(y)/7000, sd = sd(y)/1000, 
#'     stl_features(y)[c("trend", "seasonal_strength")])
#' }
#' y <- simulate_target(
#'   length = length(USAccDeaths),
#'   seasonal_periods = frequency(USAccDeaths),
#'   feature_function = my_features, target = my_features(USAccDeaths)
#' )
#' tsp(y) <- tsp(USAccDeaths)
#' plot(cbind(USAccDeaths, y))
#' cbind(my_features(USAccDeaths), my_features(y))
#' }
#' @export
#'
simulate_target <- function(length, seasonal_periods = 1, feature_function, target,
    p = 2, k = ifelse(length(seasonal_periods)==1, 3, length(seasonal_periods)),
    tolerance = 0.1, parallel = FALSE) {
  # Test lengths
  feature_length <- length(feature_function(forecast::msts(rnorm(length), 
      seasonal.periods = unique(seasonal_periods))))
  if (length(target) != feature_length) {
    stop(paste("target should be of length", feature_length))
  }

  # How many components to use? 
  if (length(seasonal_periods) == 1L) {
    k <- 3
  } else {
    k <- length(seasonal_periods)
  }
  # Set min and max for all parameters
  ga_min <- rep(0, (p + 3) * k - 1)
  ga_max <- rep(1, (p + 3) * k - 1)
  # AR parameters
  ga_min[seq((p + 1) * k)] <- -min(3,p)
  ga_max[seq((p + 1) * k)] <- min(3,p)
  # phi0
  ga_min[seq(k)*(p+1) - p] <- -10
  ga_max[seq(k)*(p+1) - p] <- 10
  # phip
  ga_min[seq(k)*(p+1)] <- -1
  ga_max[seq(k)*(p+1)] <- 1
  # Sigmas
  ga_max[(p + 1) * k + seq(k)] <- 5
  
  GA <- ga_ts(
    type = "real-valued", fitness = fitness_mar,
    length = length, seasonal_periods = seasonal_periods, ncomponents = k,
    feature_function = feature_function, target = target,
    n = length,
    min = ga_min,
    max = ga_max,
    parallel = parallel, popSize = 100, maxiter = 200,
    pmutation = 0.3, pcrossover = 0.8, maxFitness = -length(target)*(tolerance^2),
    run = 100, keepBest = TRUE, monitor = GA::gaMonitor
  )
  # Final iteration from GA algorithm
  best <- as.matrix(utils::tail(GA@bestSol, 1)[[1]])
  # Find best solution
  select <- which.min(colSums(sweep(apply(best, 2, feature_function), STATS=target, MARGIN=1)^2))
  # Return as msts or ts object
  forecast::msts(best[,select], seasonal.periods = unique(seasonal_periods))
}


fitness_mar <- function(pars, length, seasonal_periods, ncomponents, feature_function, target) {
  npar <- length(pars)
  k <- ncomponents
  p <- (npar + 1) / k - 3
  ar <- matrix(pars[seq((p + 1) * k)], ncol = k)
  sigmas <- pars[(p + 1) * k + seq(k)]
  weights <- pars[(p + 2) * k + seq(k-1)]
  weights <- exp(c(weights, 1-sum(weights)))
  weights <- weights / sum(weights)
  x <- simulate(
    mar_model(ar = ar, sigmas = sigmas, weights = weights, seasonal_periods = seasonal_periods),
    nsim = length
  )
  features <- try(feature_function(x), silent=TRUE)
  if("try-error" %in% class(features) | (any(is.na(features)))) {
    value <- -Inf
    warning("try error")
  } else {
    value <- -sum((features - target)^2)
  }
  return(list(value = value, x = x))
}
