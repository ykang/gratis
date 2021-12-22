#' Generate synthetic data from a Mixture Autoregressive model
#'
#' This function simulates one random sample path from a mixture of k Gaussian AR(p) processes.
#' The model is of the form
#' \deqn{y_t = \phi_{0,i} + \phi_{1,i}y_{t-1} + \dots + \phi_{p,i}y_{t-p} + \sigma_{i,t}\epsilon_t}
#' with probability \eqn{\alpha_i}, where \eqn{\epsilon_t} is a N(0,1) variate.
#' @param object A `mar` object, usually the output of \code{\link{mar_model}()}.
#' @param nsim length of series to generate
#' @param seed Either \code{NULL} or an integer that will be used in a call to
#' \code{\link[base]{set.seed}} before simulating the time series. The default,
#' \code{NULL}, will not change the random generator state.
#' @param n.start Length of 'burn-in' period.
#' @param ... Other arguments, not currently used.
#' @return `ts` object of length \code{nsim}.
#' @references Feng Li, Mattias Villani, and Robert Kohn. (2010). Flexible Modeling of
#'     Conditional Distributions using Smooth Mixtures of Asymmetric Student T Densities,
#'     Journal of Statistical Planning and Inference, 140(12), pp. 3638-3654.
#' @author Rob J Hyndman
#' @seealso \code{\link{mar_model}}
#' @examples
#' # MAR model with constant variances
#' phi <- cbind(c(0, 0.8, 0), c(0, 0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model1 <- mar_model(phi = phi, sigmas = c(1, 2), weights = weights)
#' y <- simulate(model1, 100)
#' plot(y)
#'
#' # MAR model with heteroskedastic errors
#' sigmas.spec <- list(
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#' )
#' model2 <- mar_model(phi = phi, sigmas = sigmas.spec, weights = weights)
#' y <- simulate(model2, 100)
#' plot(y)
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
