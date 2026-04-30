#' Simulate from a mixture autoregressive model
#'
#' Simulate one random sample path from a mixture of Gaussian autoregressive
#' components created by \code{\link{mar_model}()}. At each time point, a
#' component is sampled according to the model weights, that component's
#' conditional mean and innovation variance are used, and the resulting series
#' is returned after discarding a burn-in period.
#'
#' @param object A \code{mar} object, usually created by
#'   \code{\link{mar_model}()}.
#' @param nsim Length of the generated series.
#' @param seed Either \code{NULL} or an integer that will be used in a call to
#'   \code{\link[base]{set.seed}()} before simulating the time series. The
#'   default, \code{NULL}, leaves the random number generator state unchanged.
#' @param n.start Length of the burn-in period discarded before returning the
#'   simulated series.
#' @param ... Other arguments, not currently used.
#' @return A \code{\link[forecast]{msts}} object of length \code{nsim}. For
#'   single-seasonal or non-seasonal models this also behaves like a base
#'   \code{\link[stats]{ts}} object.
#' @references Feng Li, Mattias Villani, and Robert Kohn. (2010). Flexible Modeling of
#'     Conditional Distributions using Smooth Mixtures of Asymmetric Student T Densities,
#'     Journal of Statistical Planning and Inference, 140(12), pp. 3638-3654.
#' @author Rob J Hyndman
#' @seealso \code{\link{mar_model}}
#' @examples
#' # MAR model with constant variances
#' phi <- cbind(c(0.8, 0), c(0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model1 <- mar_model(phi = phi, sigmas = c(1, 2), weights = weights)
#' y <- simulate(model1, 100)
#' plot(y)
#'
#' # MAR model with heteroskedastic errors
#' if (requireNamespace("fGarch", quietly = TRUE)) {
#'   sigmas_spec <- list(
#'     fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'     fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#'   )
#'   model2 <- mar_model(phi = phi, sigmas = sigmas_spec, weights = weights)
#'   y <- simulate(model2, 100)
#'   plot(y)
#' }
#' @export
simulate.mar <- function(object, nsim = 100, seed = NULL, n.start = 100, ...) {
  # Set seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  # Total length of simulated series
  k <- NCOL(object$ar)
  p <- NROW(object$ar) - 1L
  n <- nsim + n.start + p

  # Randomly select components for each time
  components <- sample(seq(k), size = n, replace = TRUE, prob = object$weights)

  # Convert sigmas to n x k matrix of variances
  sigmatrix <- matrix(0, ncol = k, nrow = n)
  sigmas <- numeric(n)
  if (is.list(object$sigmas)) {
    if (!requireNamespace("fGarch", quietly = TRUE)) {
      stop("Package 'fGarch' is required to simulate MAR models with GARCH specifications.", call. = FALSE)
    }
    # GARCH components
    for (i in seq(k)) {
      sigmatrix[, i] <- fGarch::garchSim(object$sigmas[[i]],
        extended = TRUE, n = n, n.start = 10
      )$sigma
    }
  } else {
    # Homoskedastic errors
    sigmatrix <- matrix(rep(object$sigmas, n), byrow = TRUE, ncol = k)
  }
  # Select variance for each observation
  for (i in seq(k)) {
    sigmas[components == i] <- sigmatrix[components == i, i]
  }
  # Generate noise
  epsilon <- rnorm(n, sd = sqrt(sigmas))

  # Generate series
  y <- rep(0, n)
  for (i in (p + 1L):n) {
    y[i] <- sum(object$ar[, components[i]] * c(1, y[i - seq(p)])) + epsilon[i]
  }
  return(forecast::msts(y[-seq(n.start + p)], seasonal.periods = unique(object$m)))
}
