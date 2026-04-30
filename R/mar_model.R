#' Specify a mixture autoregressive model
#'
#' Construct a mixture autoregressive (MAR) model made up of \eqn{k} Gaussian
#' ARIMA\eqn{(p,d,0)(P,D,0)_m} components. The resulting \code{mar} object can
#' be passed to \code{\link{simulate.mar}()} to simulate one series or to
#' \code{\link{generate.mar}()} to generate a tsibble of several series.
#'
#' Component \eqn{i} has the form
#' \deqn{
#' (1-B)^{d_i}(1-B^{m_i})^{D_i}
#' (1-\phi_i(B))(1-\Phi_i(B)) y_t =
#' c_i + \sigma_{i,t}\epsilon_t
#' }
#' and is selected with probability \eqn{\alpha_i}. Here \eqn{B} is the
#' backshift operator, \eqn{m_i} is the seasonal period, \eqn{\epsilon_t} is a
#' standard normal variate, and \eqn{\phi_i(B)} and \eqn{\Phi_i(B)} are the
#' non-seasonal and seasonal autoregressive polynomials.
#'
#' Any missing model argument is randomly selected. Randomly selected
#' non-seasonal AR orders are sampled from \eqn{\{0,1,2,3\}},
#' non-seasonal differences from \eqn{\{0,1,2\}}, seasonal AR orders from
#' \eqn{\{0,1,2\}}, and seasonal differences from \eqn{\{0,1\}}. AR parameters
#' are sampled from the stationary region, constants from \eqn{(-3,3)}, sigmas
#' from \eqn{(1,5)}, and mixture weights from \eqn{(0,1)} before normalization.
#' If \code{k} is omitted, the number of components is sampled from
#' \eqn{\{1,2,3,4,5\}}, unless it can be inferred from another argument.
#'
#' @param k Number of mixture components.
#' @param p Non-negative integer vector giving the non-seasonal AR orders
#'   \eqn{p_i}. Ignored when \code{phi} is supplied.
#' @param d Non-negative integer vector giving the non-seasonal differencing
#'   orders \eqn{d_i}.
#' @param phi Numeric matrix of non-seasonal AR parameters. Column \eqn{i}
#'   contains the parameters for component \eqn{i}.
#' @param P Non-negative integer vector giving the seasonal AR orders
#'   \eqn{P_i}. Ignored when all \code{seasonal_periods} are 1 or when
#'   \code{Phi} is supplied.
#' @param D Non-negative integer vector giving the seasonal differencing orders
#'   \eqn{D_i}. Ignored when all \code{seasonal_periods} are 1.
#' @param Phi Numeric matrix of seasonal AR parameters. Column \eqn{i}
#'   contains the parameters for component \eqn{i}. Ignored when all
#'   \code{seasonal_periods} are 1.
#' @param constants Numeric vector of length \code{k} containing the component
#'   constants \eqn{c_1,\dots,c_k}.
#' @param sigmas Numeric vector of length \code{k}, or a list of \code{k}
#'   GARCH specifications created by \code{\link[fGarch]{garchSpec}()}.
#' @param weights Positive numeric vector of length \code{k} containing the
#'   component probabilities \eqn{\alpha_1,\dots,\alpha_k}. Values are
#'   normalized to sum to 1.
#' @param seasonal_periods Scalar seasonal period applied to every component, or
#'   a numeric vector of length \code{k} giving the seasonal period for each
#'   component.
#' @return An object of class \code{mar}. The object stores the component orders,
#'   coefficients, constants, sigmas, normalized weights, seasonal periods, and
#'   the expanded AR recursion coefficients used by \code{\link{simulate.mar}()}.
#' @author Rob J Hyndman
#' @seealso \code{\link{simulate.mar}}
#' @examples
#' # Quarterly MAR model with randomly selected parameters
#' model1 <- mar_model(seasonal_periods = 4)
#'
#' # Daily MAR model with randomly selected parameters
#' model2 <- mar_model(seasonal_periods = c(7, 365))
#'
#' # MAR model with constant variances and user-specified coefficients
#' phi <- cbind(c(0.8, 0), c(0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model3 <- mar_model(phi = phi, d = 0, sigmas = c(1, 2), weights = weights)
#'
#' # MAR model with heteroskedastic errors
#' if (requireNamespace("fGarch", quietly = TRUE)) {
#'   sigmas_spec <- list(
#'     fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'     fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#'   )
#'   model4 <- mar_model(phi = phi, sigmas = sigmas_spec, weights = weights)
#' }
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
    stop("Dimension of constants does not match other components")
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
  } else if (length(seasonal_periods) != k) {
    stop("Dimension of seasonal_periods does not match other components")
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
      P <- rep(NROW(Phi), k)
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
  if (length(nonzero) > 0) {
    ar <- ar[seq(max(nonzero)), , drop = FALSE]
  } else {
    ar <- ar[1L, , drop = FALSE]
  }

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
