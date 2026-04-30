#' Generate random values from a mixture of multivariate normal distributions
#'
#' Draw random values from a finite mixture of \eqn{k} multivariate normal
#' distributions, each with dimension \eqn{q}. For univariate mixtures,
#' \code{means} and \code{sigmas} may be supplied as vectors.
#'
#' @param n Number of observations to generate.
#' @param means A \eqn{q \times k} matrix of component means, or a vector of
#'   length \eqn{k} for univariate mixtures.
#' @param sigmas A \eqn{q \times q \times k} array of component covariance
#'   matrices, or a vector of length \eqn{k} for univariate mixtures.
#' @param weights Numeric vector of component probabilities. Values are used as
#'   sampling probabilities for the mixture components.
#' @return An \eqn{n \times q} matrix of generated data, or a numeric vector
#'   when \eqn{q = 1}.
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
  # Convert vector inputs to matrix or array
  if (is.null(dim(means))) {
    means <- t(means)
    sigmas <- array(sigmas, dim = c(1, 1, length(sigmas)))
  }
  k <- length(weights) # k-components
  q <- NROW(means) # q-dimensional
  if (!identical(dim(means), c(q, k))) {
    stop("means must be a q x k matrix")
  }
  if (!identical(dim(sigmas), c(q, q, k))) {
    stop("sigmas must be a q x q x k array.")
  }

  # Random draws for which component to use for each observation
  idx <- sample(seq(k), size = n, prob = weights, replace = TRUE)
  idx <- factor(idx, levels = seq(5))
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
  # Convert vector inputs to matrix or array
  if (is.null(dim(means))) {
    means <- t(means)
    sigmas <- array(sigmas, dim = c(1, 1, length(sigmas)))
    x <- as.matrix(x, ncol = NROW(means))
  }
  k <- length(weights) # k-components
  q <- NROW(means) # q-dimensional
  if (!identical(dim(means), c(q, k))) {
    stop("means must be a q x k matrix")
  }
  if (!identical(dim(sigmas), c(q, q, k))) {
    stop("sigmas must be a q x q x k array.")
  }

  # Compute density of each component
  component_density <- apply(
    matrix(1:k), 1,
    function(comp.i, x, means, sigmas, q) {
      mvtnorm::dmvnorm(
        x = x, mean = means[, comp.i],
        sigma = as.matrix(sigmas[, , comp.i]),
        checkSymmetry = FALSE
      )
    },
    x = x, means = means, sigmas = sigmas, q = q
  )
  density <- component_density %*% matrix(weights)
  if (log == TRUE) {
    return(log(density))
  } else {
    return(density)
  }
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
#'
#' @param n Number of observations to generate.
#' @param means.ar.par.list List of component AR parameter vectors. The first
#'   element of each vector is the intercept, followed by lag coefficients.
#' @param sigmas.list List of component variance or volatility vectors, each of
#'   length \code{n}.
#' @param weights Numeric vector of component probabilities.
#' @param yinit Initial value used to seed the recursion.
#' @return A \code{\link[stats]{ts}} object of length \code{n}.
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
#'   function(x) x$sigma
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
  y <- rep(yinit, n)
  nComp <- length(means.ar.par.list)
  nLags <- lapply(means.ar.par.list, function(x) length(x) - 1)
  maxLag <- max(unlist(nLags))

  if (any(unlist(nLags) < 1)) {
    stop("Drift is always included. Set the first elements in ar.par.list be zero to remove the drift effect.")
  }

  sigmas.ary <- do.call(cbind, sigmas.list)

  for (i in (maxLag + 1):n)
  {
    meansComp <- lapply(means.ar.par.list, function(par, y, i) {
      nLag <- length(par) - 1
      yPre <- c(1, y[(i - 1):(i - nLag)])
      yCurr <- sum(par * yPre)
      return(yCurr)
    }, y = y, i = i)
    sigmas.i <- array(sigmas.ary[i, ], dim = c(1, 1, nComp))
    y[i] <- rmixnorm(
      n = 1, means = matrix(unlist(meansComp), nrow = 1),
      sigmas = sigmas.i, weights = weights
    )
  }
  return(as.ts(y))
}


dmixnorm_ts <- function(y, means.ar.par.list, sigmas.list, weights, log = FALSE) {
  nComp <- length(means.ar.par.list)
  nLags <- lapply(means.ar.par.list, function(x) length(x) - 1)
  maxLag <- max(unlist(nLags))
  n <- length(y)

  if (any(unlist(nLags) < 1)) {
    stop("Drift is always included. Set the first elements in ar.par.list be zero to remove the drift effect.")
  }

  out.log <- y
  out.log[] <- 0 # Set the density be zero for the first maxlags.

  sigmas.ary <- matrix(do.call(cbind, sigmas.list), n, nComp, byrow = TRUE)

  for (i in (maxLag + 1):n)
  {
    meansComp <- lapply(means.ar.par.list, function(par, y, i) {
      nLag <- length(par) - 1
      yPre <- c(1, y[(i - 1):(i - nLag)])
      yCurr <- sum(par * yPre)
      ## print(cbind(yPre, par))
      return(yCurr)
    }, y = y, i = i)
    sigmas.i <- array(sigmas.ary[i, ], dim = c(1, 1, nComp))
    out.log[i] <- dmixnorm(
      x = y[i], means = matrix(unlist(meansComp), nrow = 1),
      sigmas = sigmas.i, weights = weights, log = TRUE
    )
  }

  if (log == TRUE) {
    out <- out.log
  } else {
    out <- exp(log)
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
