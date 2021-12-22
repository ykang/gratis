#' Specify parameters for a Mixture Autoregressive model
#'
#' This function allows the parameters of a mixture of k Gaussian ARIMA(p,d,0)(P,D,0)[m] processes
#' to be specified. The output is used in \code{\link{simulate.mar}()} and \code{\link{generate.mar}}.
#' The model is of the form
#' \deqn{(1-B)^{d_i}(1-B^{m_i})^{D_i} (1-\phi_i(B))(1-\Phi_i(B)) y_t = c_i + \sigma_{i,t}\epsilon_t}
#' with probability \eqn{\alpha_i}, where \eqn{B} is the backshift operator,
#' \eqn{m_i} is the seasonal period, \eqn{\epsilon_t} is a N(0,1) variate, and
#' \eqn{\phi_i(B)} and \eqn{\Phi_i(B)} are polynomials in B of
#' order \eqn{d_i} and \eqn{D_i} respectively.
#' If any argument is \code{NULL}, the corresponding parameters are randomly selected.
#' When randomly selected, the AR parameters are uniformly sampled from the stationary region,
#' p is in \{0,1,2,3\}, d is in \{0,1,2\}, P is in \{0,1,2\} and D is in \{0,1\}.
#' The model orders are uniformly sampled. The constants are uniformly sampled on (-3,3).
#' The sigmas are uniformly sampled on (1,5) and the weights are uniformly sampled on (0,1).
#' The number of components is uniformly sampled on \{1,2,3,4,5\}.
#' @param k Number of components.
#' @param p Non-negative integer vector giving the orders of non-seasonal AR polynomials \eqn{\phi_i(B)}.
#' Ignored if \code{phi} provided.
#' @param d Non-negative integer vector giving the orders of non-seasonal differencing.
#' @param phi A max(p) x k numeric matrix containing the non-seasonal AR parameters
#' (\eqn{\phi_{1,i},\dots,\phi_{p,i}}), \eqn{i=1,\dots,k} for each component.
#' @param P Non-negative integer giving the orders of seasonal AR polynomiasl \eqn{\Phi_i(B)}.
#' Ignored if \code{seasonal.periods==1} or \code{Phi} provided.
#' @param D Non-negative integer giving the orders of seasonal differencing.
#' Ignored if \code{seasonal.periods==1}.
#' @param Phi A max(P) x k numeric matrix containing the seasonal AR parameters
#' (\eqn{\Phi_{1,i},\dots,\phi_{P,i}}), \eqn{i=1,\dots,k} for each component.
#' Ignored if \code{seasonal.periods==1}.
#' @param constants A numeric vector of length k containing \eqn{c_1,\dots,c_k}.
#' @param sigmas A numeric vector of length k or a list of k GARCH specifications.
#' If it is a vector, it is assumed \eqn{\sigma_{i,t} = \sigma_i} and
#' \code{sigmas} = \eqn{\sigma_1,\dots,\sigma_k}.
#' If it is a list, each element should be the output from \code{fGarch::\link[fGarch]{garchSpec}()}.
#' @param weights A numeric vector of length k containing the probability of
#' each of the component processes, \eqn{\alpha_1,\dots,\alpha_k}.
#' @param seasonal_periods Either a scalar or a numeric vector of length k containing
#' the seasonal period of each component.
#' @return A `mar` object containing a list of \code{k}, \code{m},
#' \code{p}, \code{d}, \code{P}, \code{D},
#' \code{phi}, \code{Phi}, \code{sigmas} and \code{weights}.
#' @author Rob J Hyndman
#' @seealso \code{\link{simulate.mar}}
#' @examples
#' n <- 100
#' # Quarterly MAR model with randomly selected parameters
#' model1 <- mar_model(seasonal_periods = 4)
#'
#' # Daily MAR model with randomly selected parameters
#' model2 <- mar_model(seasonal_periods = c(7, 365))
#'
#' # MAR model with constant variances
#' # containing an AR(1) component and an AR(2) component
#' phi <- cbind(c(0, 0.8, 0), c(0, 0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model3 <- mar_model(phi = phi, d = 0, sigmas = c(1, 2), weights = weights)
#'
#' # MAR model with heteroskedastic errors
#' sigmas.spec <- list(
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#' )
#' model4 <- mar_model(phi = phi, sigmas = sigmas.spec, weights = weights)
#' @export
mar_model <- function(k = NULL,
                      p = NULL, d = NULL, phi = NULL,
                      P = NULL, D = NULL, Phi = NULL,
                      constants = NULL, sigmas = NULL, weights = NULL,
                      seasonal_periods = 1L) {
  # How many components?
  if (is.null(k)) {
    if (!is.null(weights)) {
      k <- length(weights)
    } else if (!is.null(phi)) {
      k <- NCOL(phi)
    } else if (!is.null(Phi)) {
      k <- NCOL(Phi)
    } else if (!is.null(sigmas)) {
      k <- length(sigmas)
    } else if (!is.null(constants)) {
      k <- length(constants)
    } else if (length(seasonal_periods) == 1L) {
      # choose a random number of components between 1 and 5
      k <- sample(seq(5), 1)
    } else {
      k <- length(seasonal_periods)
    }
  }

  # weights
  if (is.null(weights)) {
    # Choose a random set of weights in (0,1)
    weights <- runif(k)
  } else {
    if (length(weights) != k) {
      stop("Dimension of weights does not match other components")
    } else if (any(weights <= 0)) {
      stop("weights must be positive")
    }
  }
  # Rescale weights so they sum to 1
  weights <- weights / sum(weights)

  # constants
  if (is.null(constants)) {
    # Choose k constants in [-3,3]
    constants <- runif(k, -3, 3)
  } else if (length(constants) != k) {
    stop("Dimension of sigmas does not match other components")
  }

  # sigmas
  if (is.null(sigmas)) {
    # Choose k variances in [1,5]
    sigmas <- runif(k, 1, 5)
  } else {
    if (length(sigmas) != k) {
      stop("Dimension of sigmas does not match other components")
    } else if (is.numeric(sigmas)) {
      if (any(sigmas <= 0)) {
        stop("sigmas must be positive")
      }
    }
  }
  # one seasonal period for each component
  if (length(seasonal_periods) == 1L) {
    seasonal_periods <- rep(seasonal_periods, k)
  }
  # Non-seasonal AR parameters
  if (is.null(p) & is.null(phi)) {
    p <- sample(c(0, 1, 2, 3), k, replace = TRUE)
  } else if (!is.null(phi)) {
    if (NCOL(phi) != k) {
      stop("Dimension of phi does not match other components")
    }
    p <- rep(NROW(phi), k)
  }
  if (is.null(phi)) {
    if (length(p) == 1L) {
      p <- rep(p, k)
    }
    phi <- matrix(0, ncol = k, nrow = max(p))
    for (i in seq(k)) {
      if (p[i] > 0) {
        phi[seq(p[i]), i] <- stationary_ar(p[i])
      }
    }
  }
  if (is.null(d)) {
    d <- sample(c(0, 1, 2), k, replace = TRUE)
  } else if (length(d) == 1L) {
    d <- rep(d, k)
  }
  # Seasonal AR parameters
  if (any(seasonal_periods > 1)) {
    if (is.null(P) & is.null(Phi)) {
      P <- sample(c(0,1,2), k, replace=TRUE)
    } else if (!is.null(Phi)) {
      if (NCOL(Phi) != k) {
        stop("Dimension of Phi does not match other components")
      }
      P <- rep(NROW(phi), k)
    }
    if (is.null(Phi)) {
      if (length(P) == 1L) {
        P <- rep(P, k)
      }
      P[seasonal_periods == 1] <- 0
      Phi <- matrix(0, ncol = k, nrow = max(P))
      for (i in seq(k)) {
        if (P[i] > 0) {
          Phi[seq(P[i]), i] <- stationary_ar(P[i])
        }
      }
    }
    if (is.null(D)) {
      D <- sample(c(0, 1), k, replace = TRUE)
    } else if (length(D) == 1L) {
      D <- rep(D, k)
    }
    D[seasonal_periods == 1 | d == 2] <- 0
  } else {
    D <- P <- rep(0, k)
    Phi <- matrix(0, ncol = k, nrow = 0)
  }
  constants[d + D > 1] <- 0
  # Compute pi weights
  max_order <- max(d + seasonal_periods * D + p + seasonal_periods * P)
  ar <- matrix(0, ncol = k, nrow = max_order + 1)
  for (i in seq(k)) {
    pi_weights <- c(
      constants[i],
      pi_coefficients(ar = phi[, i], sar = Phi[, i], d = d[i], D = D[i], m = seasonal_periods[i])
    )
    ar[seq_along(pi_weights), i] <- pi_weights
  }
  # Trim redundant zeros
  nonzero <- which(rowSums(abs(ar)) > 0)
  ar <- ar[seq(max(nonzero)), , drop = FALSE]

  # Return mar object
  structure(list(
    k = k, p = p, d = d, phi = phi, P = P, D = D, Phi = Phi,
    constants = constants, ar = ar,
    sigmas = sigmas, weights = weights, m = seasonal_periods
  ),
  class = "mar"
  )
}

#' @method print mar
#' @export
print.mar <- function(x, ...) {
  cat(paste("Mixture AR model with",x$k,"components:\n"))
  for(i in seq(x$k)) {
    cat(paste0("    ARIMA(",x$p[i],",",x$d[i],",","0)"))
    if(x$P[i] > 0 | x$D[i] > 0) {
      cat(paste0("(",x$P[i],",",x$D[i],",","0)[",x$m[i],"]"))
    }
    cat(" with weight ")
    cat(sprintf("%.2f",round(x$weights[i],2)))
    cat("\n")
  }
}
