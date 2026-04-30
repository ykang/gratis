#' Generate time series with target features
#'
#' Search for MAR model parameters that produce time series with feature values
#' close to a user-specified target. \code{simulate_target()} returns one
#' generated series, while \code{generate_target()} repeats the search to return
#' several generated series in a tsibble.
#'
#' A genetic algorithm evaluates candidate MAR parameter vectors by simulating a
#' series, scaling it, computing \code{feature_function()}, and comparing those
#' features to \code{target}. Because each candidate series is scaled before
#' feature evaluation, the target features should not depend on the absolute
#' scale of the series.
#'
#' @param length Length of each generated time series.
#' @param seasonal_periods Scalar seasonal period, or a numeric vector of
#'   seasonal periods for multiple-seasonal data.
#' @param feature_function Function that accepts a time series and returns a
#'   numeric vector of features.
#' @param target Numeric vector of target feature values. It must have the same
#'   length and order as the vector returned by \code{feature_function()}.
#' @param k Number of MAR components to use. The default is 3 for single-seasonal
#'   or non-seasonal data, and the number of seasonal periods for
#'   multiple-seasonal data.
#' @param tolerance Average absolute feature tolerance. Larger values usually
#'   return faster but less precise results.
#' @param trace Logical value indicating whether to show genetic algorithm
#'   progress.
#' @param parallel Logical value or parallel backend specification passed to the
#'   genetic algorithm. Use \code{TRUE} to run fitness evaluations in parallel.
#' @return \code{simulate_target()} returns a \code{\link[forecast]{msts}}
#'   object of length \code{length}. \code{generate_target()} returns a tsibble
#'   containing \code{nseries} generated series.
#' @author Yanfei Kang and Rob J Hyndman
#' @examples
#' \dontrun{
#' set.seed(1)
#' library(tsfeatures)
#' my_features <- function(y) {
#'   c(entropy(y), acf = acf(y, plot = FALSE)$acf[2:3, 1, 1])
#' }
#' # Simulate a ts with specified target features
#' y <- simulate_target(
#'   length = 60, feature_function = my_features, target = c(0.5, 0.9, 0.8)
#' )
#' my_features(y)
#' plot(y)
#'
#' # Generate a tsibble with specified target features
#' df <- generate_target(
#'   length = 60, feature_function = my_features, target = c(0.5, 0.9, 0.8)
#' )
#' df %>% 
#'  as_tibble() %>%
#'  group_by(key) %>%
#'  dplyr::reframe(
#'            value = my_features(value), 
#'            feature=c("entropy","acf1", "acf2")
#'            )
#' autoplot(df)
#' # Simulate time series similar to an existing series
#' my_features <- function(y) {
#'   c(stl_features(y)[c("trend", "seasonal_strength", "peak", "trough")])
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
simulate_target <- function(length=100, seasonal_periods = 1, feature_function, target,
                            k = ifelse(length(seasonal_periods) == 1, 3, length(seasonal_periods)),
                            tolerance = 0.05, trace = FALSE, parallel = FALSE) {
  # How many components to use?
  k <- ifelse(length(seasonal_periods) == 1L, 3, length(seasonal_periods))
  # Set AR orders
  p <- 3
  P <- ifelse(max(seasonal_periods) > 1, 2, 0)
  # Parameter vector: d, D, phi, Phi, constants, sigmas, weights
  par_length <- (2+p+P+3)*k - 1
  # Set min and max for all parameters
  ga_min <- rep(0, par_length)
  ga_max <- rep(1, par_length)
  # d
  ga_max[seq(k)] <- 2
  # phi and Phi
  ga_min[seq(p*k + P*k) + 2*k] <- -1.5
  ga_max[seq(p*k + P*k) + 2*k] <- 1.5
  # constants
  ga_min[seq(k) + (p+P+2)*k] <- -3
  ga_max[seq(k) + (p+P+2)*k] <- 3
  # sigmas
  ga_max[seq(k) + (p+P+3)*k] <- 5

  if(trace)
    monitor <- GA::gaMonitor
  else
    monitor <- NULL
  GA <- ga_ts(
    type = "real-valued", fitness = fitness_mar,
    length = length, seasonal_periods = seasonal_periods, ncomponents = k,
    feature_function = feature_function, target = target,
    n = length,
    min = ga_min,
    max = ga_max,
    parallel = parallel, popSize = 100, maxiter = 200,
    pmutation = 0.1, pcrossover = 0.8, maxFitness = -tolerance,
    run = 50, keepBest = TRUE, monitor = monitor
  )
  # Final iteration from GA algorithm
  best <- forecast::msts(utils::tail(GA@bestSol, 1)[[1]],
    seasonal.periods = unique(seasonal_periods)
  )
  # Find best solution
  if (NCOL(best) > 1) {
    best_features <- matrix(0, ncol=NCOL(best), nrow=length(target))
    for(i in seq(NCOL(best))) {
      best_features[,i] <- feature_function(best[,i])
    }
    select <- which.min(colSums(sweep(best_features, STATS = target, MARGIN = 1)^2))
    best <- best[, select]
  }
  # Return as msts or ts object
  best
}

#' @rdname simulate_target
#' @param nseries Number of series to generate
#' @export
#'
generate_target <-  function(length=100, nseries=10, seasonal_periods = 1, feature_function, target,
                             k = ifelse(length(seasonal_periods) == 1, 3, length(seasonal_periods)),
                             tolerance = 0.05, trace = FALSE, parallel = FALSE) {

  tsmatrix <- matrix(NA_real_, nrow=length, ncol=nseries)
  for(i in seq(nseries)) {
    tsmatrix[,i] <- simulate_target(
      length=length, 
      seasonal_periods = seasonal_periods, 
      feature_function = feature_function, 
      target = target,
      k=k, 
      tolerance=tolerance, 
      trace=trace, 
      parallel=parallel
    )
  }
  make_tsibble(tsmatrix, seasonal_periods)
}

fitness_mar <- function(pars, length, seasonal_periods, ncomponents, feature_function, target) {
  npar <- length(pars)
  k <- ncomponents
  # AR orders
  p <- 3
  P <- ifelse(max(seasonal_periods) > 1, 2, 0)
  # Parameter vector: d, D, phi, Phi, constants, sigmas, weights
  d <- round(pars[seq(k)])
  D <- round(pars[k+seq(k)])
  phi <- matrix(pars[2*k + seq(p*k)], ncol=k, nrow=p)
  if(P > 0)
    Phi <- matrix(pars[(2+p)*k + seq(P*k)], ncol=k, nrow=P)
  else
    Phi <- matrix(0, ncol=k, nrow=0)
  constants <- pars[(2+p+P)*k + seq(k)]
  sigmas <- pars[(3+p+P)*k + seq(k)]
  weights <- pars[(4+p+P)*k + seq(k-1)]
  weights <- exp(c(weights, 1 - sum(weights)))
  weights <- weights / sum(weights)
  x <- simulate(
    mar_model(d=d, D=D, p=p, P=P, phi=phi, Phi=Phi, constants=constants, 
              sigmas = sigmas, weights = weights, 
              seasonal_periods = seasonal_periods),
    nsim = length
  )
  x <- scalets(x)
  features <- try(feature_function(x), silent = TRUE)
  if ("try-error" %in% class(features) | (any(is.na(features)))) {
    value <- -Inf
    warning("try error")
  } else {
    value <- -mean(abs(features - target))
  }
  return(list(value = value, x = x))
}
