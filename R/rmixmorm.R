#' Generate random values from a mixture of multivariate normal distributions
#'
#' Draw random values from a finite mixture of \eqn{k} multivariate normal
#' distributions, each with dimension \eqn{q}. For univariate mixtures,
#' \code{means} and \code{sigmas} may be supplied as vectors.
#' Covariance matrices must be finite, symmetric and positive definite,
#' including those of zero-weight components.
#'
#' @param n Non-negative integer number of observations to generate.
#' @param means A \eqn{q \times k} matrix of component means, or a vector of
#'   length \eqn{k} for univariate mixtures.
#' @param sigmas A \eqn{q \times q \times k} array of component covariance
#'   matrices, or a vector of length \eqn{k} for univariate mixtures. These are
#'   variances or covariances, not standard deviations.
#' @param weights Finite, non-negative component weights, with at least one
#'   positive value. Weights are normalized to sum to one.
#' @return An \eqn{n \times q} matrix of generated data, or a numeric vector
#'   when \eqn{q = 1}.
#' @seealso \code{\link{rmixnorm_ts}} for autoregressive conditional means.
#' @references Villani et al 2009.
#' @author Feng Li, Central University of Finance and Economics.
#' @examples
#' out <- rmixnorm(
#'   n = 1000, means = c(-5, 0, 5), sigmas = c(1, 1, 3),
#'   weights = c(0.3, 0.4, 0.3)
#' )
#' hist(out, breaks = 100, freq = FALSE)
#' @export
rmixnorm <- function(n, means, sigmas, weights) {
  if (length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
    stop("n must be a non-negative integer")
  }
  inputs <- prepare_mixture(means, sigmas, weights)
  means <- inputs$means
  sigmas <- inputs$sigmas
  weights <- inputs$weights
  k <- NCOL(means)
  q <- NROW(means)

  # Random draws for which component to use for each observation
  idx <- sample(seq(k), size = n, prob = weights, replace = TRUE)
  idx <- factor(idx, levels = seq(k))
  nsamp <- table(idx)
  data <- matrix(NA_real_, nrow = n, ncol = q)
  # Generate draws from each component
  for (i in seq(k)) {
    if (nsamp[i] > 0) {
      data[idx == i, ] <- mvtnorm::rmvnorm(
        n = nsamp[i], mean = means[, i],
        sigma = as.matrix(sigmas[, , i]),
        checkSymmetry = FALSE
      )
    }
  }
  # Return a vector if q=1, or a matrix otherwise
  if (q == 1) {
    data <- c(data)
  }
  data
}

# Density of mixture of multivariate normals
dmixnorm <- function(x, means, sigmas, weights, log = FALSE) {
  inputs <- prepare_mixture(means, sigmas, weights)
  means <- inputs$means
  sigmas <- inputs$sigmas
  weights <- inputs$weights
  q <- NROW(means)
  x <- if (q == 1L) matrix(x, ncol = 1L) else if (is.null(dim(x))) matrix(x, nrow = 1L) else x
  if (NCOL(x) != q) stop("x must have q columns")
  components <- vapply(seq_along(weights), function(j) {
    if (weights[j] == 0) return(rep(-Inf, NROW(x)))
    log(weights[j]) + mvtnorm::dmvnorm(
      x, mean = means[, j], sigma = matrix(sigmas[, , j], q, q), log = TRUE)
  }, numeric(NROW(x)))
  components <- matrix(components, nrow = NROW(x), ncol = length(weights))
  density <- matrix(row_logsumexp(components), ncol = 1L)
  if (log) density else exp(density)
}

normalise_mixture_weights <- function(weights, k) {
  if (length(weights) != k || any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("weights must be finite and non-negative, with at least one positive")
  }
  weights <- weights / max(weights)
  weights / sum(weights)
}

prepare_mixture <- function(means, sigmas, weights) {
  if (is.null(dim(means))) {
    means <- t(means)
    sigmas <- array(sigmas, dim = c(1, 1, length(sigmas)))
  }
  k <- length(weights)
  q <- NROW(means)
  if (!identical(dim(means), c(q, k))) {
    stop("means must be a q x k matrix")
  }
  if (!identical(dim(sigmas), c(q, q, k))) {
    stop("sigmas must be a q x q x k array.")
  }

  if (q < 1L || k < 1L || any(!is.finite(means))) stop("means must be nonempty and finite")
  for (j in seq_len(k)) {
    covariance <- matrix(sigmas[, , j], q, q)
    symmetric <- all(abs(covariance - t(covariance)) <= 1e-12 * abs(t(covariance)))
    if (any(!is.finite(covariance)) || !isTRUE(symmetric) ||
        inherits(try(chol(covariance), silent = TRUE), "try-error")) {
      stop("sigma matrices must be finite, symmetric and positive definite")
    }
  }
  list(means = means, sigmas = sigmas, weights = normalise_mixture_weights(weights, k))
}

row_logsumexp <- function(values) {
  if (NROW(values) == 0L) return(numeric(0))
  largest <- values[, 1L]
  for (j in seq_len(NCOL(values))[-1L]) largest <- pmax(largest, values[, j])
  result <- largest
  finite <- is.finite(largest)
  result[finite] <- largest[finite] + log(rowSums(exp(values[finite, , drop = FALSE] - largest[finite])))
  result
}

prepare_ts_mixture <- function(n, parameters, sigmas, weights) {
  k <- length(parameters)
  if (!k || length(sigmas) != k) stop("sigmas.list must match a nonempty means.ar.par.list")
  if (any(lengths(parameters) < 2L)) stop("each AR vector must include an intercept and at least one lag")
  if (any(!is.finite(unlist(parameters)))) stop("AR parameters must be finite")
  if (any(lengths(sigmas) != n)) stop("each sigma vector must have length n")
  variances <- matrix(unlist(sigmas), nrow = n, ncol = k)
  if (any(!is.finite(variances)) || any(variances <= 0)) stop("sigmas must be positive and finite")
  list(scales = sqrt(variances), lags = lapply(lengths(parameters) - 1L, seq_len),
       weights = normalise_mixture_weights(weights, k))
}

## Tests
## n <- 1000
## means = c(-5,0,5)
## sigmas = c(1,1,3)
## weights=c(0.3,0.4,0.3)
## out <- rmixnorm(n, means, sigmas, weights)
## out.density <- dmixnorm(out, means, sigmas, weights)
## hist(out, breaks = 100, freq = FALSE)
## points(out, out.density)

#' Simulate a mixture-normal autoregressive process
#'
#' Simulate a univariate process whose conditional distribution is a finite
#' mixture of Gaussian distributions. The conditional mean of each component is
#' an autoregressive process with its own parameter vector.
#' All components use the common observed past, rather than independent latent
#' component paths. No burn-in is removed: the first \code{min(n, max_order)}
#' observations are fixed at \code{yinit}, where \code{max_order} is the largest
#' component AR order. A shorter series consists entirely of these initial values.
#'
#' @param n Non-negative integer total path length, including initial observations.
#' @param means.ar.par.list List of component AR parameter vectors. The first
#'   element of each vector is the intercept, followed by at least one lag
#'   coefficient. Orders may differ between components; values must be finite.
#' @param sigmas.list List of positive, finite component variance vectors, each of
#'   length \code{n}. Values are variances, not standard deviations.
#' @param weights Finite, non-negative component weights. At least one must be
#'   positive; weights are normalized to sum to one.
#' @param yinit Finite scalar assigned to the initial observations.
#' @return A \code{\link[stats]{ts}} object of length \code{n}, or
#'   \code{numeric(0)} for zero length. Stationarity is not enforced.
#' @seealso \code{\link{rmixnorm}}, \code{\link{mar_model}}, \code{\link{simulate.mar}}
#' @references Feng Li, Mattias Villani, and Robert Kohn. (2010). Flexible Modeling of
#'     Conditional Distributions using Smooth Mixtures of Asymmetric Student T Densities,
#'     Journal of Statistical Planning and Inference, 140(12), pp. 3638-3654.
#' @author Feng Li, Central University of Finance and Economics.
#' @examples
#' if (requireNamespace("fGarch", quietly = TRUE)) {
#' n <- 1000
#' means.ar.par.list <- list(c(0, 0.8), c(0, 0.6, 0.3))
#' sigmas.spec <- list(
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.06)), cond.dist = "norm"),
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)), cond.dist = "norm")
#' )
#' sigmas.list <- lapply(
#'   lapply(sigmas.spec, fGarch::garchSim, extended = TRUE, n = n),
#'   function(x) x$sigma^2
#' )
#' weights <- c(0.8, 0.2)
#' y <- rmixnorm_ts(
#'   n = n, means.ar.par.list = means.ar.par.list, sigmas.list = sigmas.list,
#'   weights = weights
#' )
#' plot(y)
#' }
#' @export
rmixnorm_ts <- function(n, means.ar.par.list, sigmas.list, weights, yinit = 0) {
  if (length(n) != 1L || !is.finite(n) || n < 0 || n != as.integer(n)) {
    stop("n must be a non-negative integer")
  }
  inputs <- prepare_ts_mixture(n, means.ar.par.list, sigmas.list, weights)
  if (length(yinit) != 1L || !is.finite(yinit)) stop("yinit must be finite")
  y <- rep(yinit, n)
  nComp <- length(means.ar.par.list)
  maxLag <- max(lengths(inputs$lags))

  for (i in seq_len(max(0L, n - maxLag)) + maxLag)
  {
    component <- sample.int(nComp, 1L, prob = inputs$weights)
    mean <- sum(means.ar.par.list[[component]] * c(1, y[i - inputs$lags[[component]]]))
    y[i] <- rnorm(1L, mean, inputs$scales[i, component])
  }
  if (n == 0L) return(numeric(0))
  return(as.ts(y))
}


dmixnorm_ts <- function(y, means.ar.par.list, sigmas.list, weights, log = FALSE) {
  nComp <- length(means.ar.par.list)
  n <- length(y)
  inputs <- prepare_ts_mixture(n, means.ar.par.list, sigmas.list, weights)
  maxLag <- max(lengths(inputs$lags))

  out.log <- y
  out.log[] <- 0 # Set the density be zero for the first maxlags.

  if (n > maxLag) {
    rows <- seq_len(n - maxLag) + maxLag
    values <- matrix(-Inf, nrow = length(rows), ncol = nComp)
    past <- c(0, utils::head(as.numeric(y), -1L))
    # Conditional means depend only on observed lags, so each component can
    # be evaluated for the entire path with one finite filter.
    for (j in which(inputs$weights > 0)) {
      parameters <- means.ar.par.list[[j]]
      means <- parameters[1L] + stats::filter(past, parameters[-1L], sides = 1)
      values[, j] <- log(inputs$weights[j]) +
        stats::dnorm(y[rows], means[rows], inputs$scales[rows, j], log = TRUE)
    }
    out.log[rows] <- row_logsumexp(values)
  }

  if (log == TRUE) {
    out <- out.log
  } else {
    out <- exp(out.log)
  }
  return(out)
}

## n = 1000
## means.ar.par.list = list(c(0, 0.8), c(0, 0.6, 0.3))

## require("fGarch")
## sigmas.spec <- list(fGarch::garchSpec(model = list(alpha = c(0.05, 0.06)), cond.dist = "norm"),
##                     fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)), cond.dist = "norm"))
## sigmas.list <- lapply(lapply(sigmas.spec, fGarch::garchSim, extended = TRUE, n = n),
##                       function(x) x$sigma)
## weights <- c(0.8, 0.2)
## y = rmixnorm_ts(n = n, means.ar.par.list = means.ar.par.list, sigmas.list = sigmas.list,
##                 weights = weights)
## out = dmixnorm_ts(y = y, means.ar.par.list = means.ar.par.list,
##                   sigmas.list = sigmas.list, weights = weights, log = TRUE)
