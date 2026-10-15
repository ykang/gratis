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
#' All components condition on the same generated past; the result is not a
#' weighted average of independently simulated component paths. Numeric
#' \code{object$sigmas} are innovation variances, not standard deviations.
#'
#' @param nsim Non-negative integer output length, excluding burn-in.
#' @param seed Optional integer passed to \code{\link[base]{set.seed}()} before
#'   simulation. \code{NULL} uses and advances the current random stream without
#'   resetting it. Use the \code{seed} argument to \code{\link{generate.mar}()}
#'   to seed a batch once rather than restarting the stream for each path.
#' @param n.start Non-negative integer burn-in length, discarded before returning the
#'   simulated series. The default of 100 does not guarantee convergence to a
#'   stationary distribution, particularly for near-unit-root models. Stationary
#'   components do not by themselves guarantee mixture stability.
#' @param ... Other arguments, not currently used.
#' @return A \code{\link[forecast]{msts}} object of length \code{nsim}. For
#'   single-seasonal or non-seasonal models this also behaves like a base
#'   \code{\link[stats]{ts}} object. Zero output length returns \code{numeric(0)}.
#'   Values are not centered or standardized.
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
  for (value in list(nsim, n.start)) {
    if (length(value) != 1L || !is.finite(value) || value < 0 || value != as.integer(value)) {
      stop("nsim and n.start must be non-negative integers")
    }
  }
  if (nsim == 0L) return(numeric(0))
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

  if (is.list(object$sigmas)) {
    sigmatrix <- matrix(0, ncol = k, nrow = n)
    if (!requireNamespace("fGarch", quietly = TRUE)) {
      stop("Package 'fGarch' is required to simulate MAR models with GARCH specifications.", call. = FALSE)
    }
    # GARCH components
    for (i in seq(k)) {
      sigmatrix[, i] <- fGarch::garchSim(object$sigmas[[i]],
        extended = TRUE, n = n, n.start = 10
      )$sigma^2
    }
    sigmas <- sigmatrix[cbind(seq_len(n), components)]
  } else {
    # Homoskedastic errors
    sigmas <- object$sigmas[components]
  }
  # Generate noise
  epsilon <- rnorm(n, sd = sqrt(sigmas))

  # Generate series
  if (p == 0L) {
    y <- object$ar[1L, components] + epsilon
  } else if (k == 1L) {
    y <- stats::filter(c(rep(0, p), object$ar[1L, 1L] + epsilon[seq.int(p + 1L, n)]),
                       object$ar[-1L, 1L], method = "recursive")
  } else {
    y <- rep(0, n)
    lags <- lapply(seq_len(k), function(j) which(object$ar[-1L, j] != 0))
    coefficients <- lapply(seq_len(k), function(j) object$ar[lags[[j]] + 1L, j])
    for (i in (p + 1L):n) {
      component <- components[i]
      y[i] <- object$ar[1L, component] +
        sum(coefficients[[component]] * y[i - lags[[component]]]) + epsilon[i]
    }
  }
  return(forecast::msts(as.numeric(utils::tail(y, nsim)), seasonal.periods = unique(object$m)))
}
