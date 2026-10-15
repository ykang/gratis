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
#' Finite nonconstant series are centered and sample-standardized using
#' overflow-safe scaling; exactly constant series are left unchanged.
#'
#' Mixture weights are the spacings between sorted breakpoints in \eqn{[0,1]},
#' allowing the full simplex, including zero-weight components. If tolerance
#' is not reached, the best series is returned with a warning giving its mean
#' absolute feature error. Callback errors and incompatible feature lengths
#' are reported as errors; nonfinite candidate paths or features are rejected.
#' The best final population member is selected using its cached fitness and
#' series; generation histories are not retained or callbacks reevaluated.
#'
#' @param length Positive integer output length of each generated series.
#' @param seasonal_periods Positive integer seasonal period, applied to every
#'   component, or a vector of length \code{k} specifying one period per component.
#' @param feature_function Function that accepts a candidate time series and
#'   returns one numeric value per target feature. It should be deterministic
#'   for reproducible scores. Nonfinite values reject a candidate; errors raised
#'   by the function propagate. Callbacks used in parallel must be available
#'   to the workers along with their dependencies.
#' @param target Numeric vector of target feature values. It must have the same
#'   length and order as the vector returned by \code{feature_function()}.
#' @param k Number of MAR components to use. The default is 3 for single-seasonal
#'   or non-seasonal data, and the number of seasonal periods for
#'   multiple-seasonal data.
#' @param tolerance Non-negative finite mean absolute feature-error tolerance.
#'   Reaching this value stops the search. Success is not guaranteed; a miss
#'   returns the best series with a warning reporting its error.
#' @param trace Logical value indicating whether to show genetic algorithm
#'   progress.
#' @param parallel Logical value or parallel backend specification passed to the
#'   genetic algorithm. Use \code{TRUE} to run fitness evaluations in parallel.
#'   A supplied cluster requires \pkg{doParallel} and remains caller-owned;
#'   stop it explicitly when it is no longer needed. Internally created clusters
#'   are stopped automatically.
#' @return \code{simulate_target()} returns a \code{\link[forecast]{msts}}
#'   object of length \code{length}. \code{generate_target()} returns a tsibble
#'   containing \code{nseries} generated series. Returned values are on the
#'   standardized scale used for feature evaluation, except that exactly constant
#'   candidates are unchanged. Neither function returns fitted MAR parameters.
#'   A search with no finite valid candidate stops with an error.
#' @seealso \code{\link{mar_model}}, \code{\link{generate.mar}}
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
  if (is.null(k)) {
    k <- ifelse(length(seasonal_periods) == 1L, 3, length(seasonal_periods))
  }
  if (!is.function(feature_function)) stop("feature_function must be a function")
  if (!is.numeric(seasonal_periods) || !length(seasonal_periods) ||
      any(!is.finite(seasonal_periods)) || any(seasonal_periods < 1) ||
      any(seasonal_periods != as.integer(seasonal_periods))) {
    stop("seasonal_periods must contain positive integers")
  }
  if (!is.numeric(target) || !length(target) || any(!is.finite(target))) {
    stop("target must be a nonempty finite numeric vector")
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be non-negative and finite")
  }
  for (value in list(length, k)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != as.integer(value)) {
      stop("length and k must be positive integers")
    }
  }
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
    run = 50, keepBest = FALSE, monitor = monitor
  )
  candidates <- GA@popTs
  scores <- GA@fitness
  scores[colSums(!is.finite(candidates)) > 0L] <- -Inf
  if (!length(scores) || !any(is.finite(scores))) {
    stop("target search did not produce a finite series")
  }
  select <- which.max(scores)
  best <- forecast::msts(candidates[, select], seasonal.periods = unique(seasonal_periods))
  error <- -scores[select]
  if (error > tolerance) {
    warning(sprintf("target tolerance not reached: mean absolute feature error=%.6g, tolerance=%.6g",
                    error, tolerance), call. = FALSE)
  }
  best
}

#' @rdname simulate_target
#' @param nseries Positive integer number of separate target searches. Searches
#'   are repeated sequentially, with successive draws from the current RNG
#'   stream; \code{parallel} controls candidate evaluations within each search.
#' @export
#'
generate_target <-  function(length=100, nseries=10, seasonal_periods = 1, feature_function, target,
                             k = ifelse(length(seasonal_periods) == 1, 3, length(seasonal_periods)),
                             tolerance = 0.05, trace = FALSE, parallel = FALSE) {
  for (value in list(length, nseries)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      stop("length and nseries must be positive integers")
    }
  }
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
  sigmas <- pmax(pars[(3+p+P)*k + seq_len(k)], 1e-8)
  weights <- diff(sort(c(0, pars[(4+p+P)*k + seq_len(k-1L)], 1)))
  x <- simulate(
    mar_model(d=d, D=D, p=p, P=P, phi=phi, Phi=Phi, constants=constants, 
              sigmas = sigmas, weights = weights, 
              seasonal_periods = seasonal_periods),
    nsim = length
  )
  if (any(!is.finite(x))) {
    return(list(value = -Inf, x = rep(NA_real_, length)))
  }
  x <- scalets(x)
  list(value = target_score(x, feature_function, target), x = x)
}

target_score <- function(x, feature_function, target) {
  if (any(!is.finite(x))) return(-Inf)
  features <- feature_function(x)
  if (!is.numeric(features) || length(features) != length(target)) {
    stop("feature_function must return one numeric value per target feature")
  }
  if (any(!is.finite(features))) return(-Inf)
  -mean(abs(features - target))
}
