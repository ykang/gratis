#' Specify parameters for a Mixture Autoregressive model
#'
#' This function allows the parameters of a mixture of k Gaussian AR processes
#' to be specified. The output is used in \code{\link{simulate.mar}()}.
#' The model is of the form
#' \deqn{y_t = \phi_{0,i} + \phi_{1,i}y_{t-1} + \cdots + \phi_{p,i}y_{t-p} + \sigma_{i,t}\epsilon_t}
#' with probability \eqn{\alpha_i}, where \eqn{\epsilon_t} is a N(0,1) variate.
#' If any argument is \code{NULL}, the corresponding parameters are randomly selected.
#' The AR parameters may be non-stationary. When randomly selected, they are chosen
#' to be the pi weights of an ARIMA(p,d,0)(P,D,0) process where p is in \{0,1,2,3\},
#' d is in \{0,1,2\}, P is in \{0,1,2\} and D is in \{0,1\}. The model orders and the
#' parameters are uniformly sampled. The sigmas are uniformly sampled on (1,5)
#' and the weights are uniformly sampled on (0,1). The number of components is
#' uniformly sampled on \{1,2,3,4,5\}.
#' @param ar A (p+1) x k numeric matrix containing the AR parameters
#' (\eqn{\phi_{0,i}, \phi_{1,i},\dots,\phi_{p,i}}), \eqn{i=1,\dots,k} for each component.
#' @param sigmas A numeric vector of length k or a list of k GARCH specifications.
#' If it is a vector, it is assumed \eqn{\sigma_{i,t} = \sigma_i} and
#' \code{sigmas} = \eqn{\sigma_1,\dots,\sigma_k}.
#' If it is a list, each element should be the output from \code{fGarch::\link[fGarch]{garchSpec}()}.
#' @param weights A numeric vector of length k containing the probability of
#' each of the component processes, \eqn{\alpha_1,\dots,\alpha_k}.
#' @param seasonal_periods Either a scalar or a numeric vector of length k containing
#' the number of seasonal periods for each component.
#' @return A `mar` object containing a list of \code{k}, \code{p}, \code{ar},
#' \code{sigmas} and \code{weights}.
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
#' ar <- cbind(c(0, 0.8, 0), c(0, 0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model3 <- mar_model(ar = ar, sigmas = c(1, 2), weights = weights)
#'
#' # MAR model with heteroskedastic errors
#' sigmas.spec <- list(
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.06))),
#'   fGarch::garchSpec(model = list(alpha = c(0.05, 0.05)))
#' )
#' model4 <- mar_model(ar = ar, sigmas = sigmas.spec, weights = weights)
#' @export
mar_model <- function(ar = NULL, sigmas = NULL, weights = NULL, seasonal_periods = 1L) {
  # How many components?
  if (is.null(ar) & is.null(sigmas) & is.null(weights) & length(seasonal_periods) == 1L) {
    # choose a random number of components between 1 and 5
    k <- sample(seq(5), 1, replace = TRUE)
  } else if (!is.null(weights)) {
    k <- length(weights)
  } else if (!is.null(sigmas)) {
    k <- length(sigmas)
  } else if (length(seasonal_periods) > 1L) {
    k <- length(seasonal_periods)
  } else {
    k <- NCOL(ar)
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
  # sigmas
  if (is.null(sigmas)) {
    # Choose k variances in [1,5]
    sigmas <- runif(k, 1, 5)
  } else {
    if (length(sigmas) != k) {
      stop("Dimension of sigmas does not match other components")
    } else if(is.numeric(sigmas)) {
      if(any(sigmas <= 0)) {
        stop("sigmas must be positive")
      }
    }
  }
  # one seasonal period for each component
  if (length(seasonal_periods) == 1L) {
    seasonal_periods <- rep(seasonal_periods, k)
  }
  # AR parameters
  if (!is.null(ar)) {
    if (NCOL(ar) != k) {
      stop("Dimension of ar does not match other components")
    }
  } else {
    ar <- matrix(0, ncol = k, nrow = 4 * max(seasonal_periods) + 4)
    p <- sample(c(0, 1, 2, 3), k, replace = TRUE)
    for (i in seq(k)) {
      if (seasonal_periods[i] > 1) {
        d <- sample(c(0, 1), 1, replace = TRUE)
        D <- sample(c(0, 1), 1, replace = TRUE)
        P <- sample(c(0, 1, 2), 1, replace = TRUE)
      } else {
        D <- P <-
          d <- sample(c(0, 1, 2), 1, replace = TRUE)
      }
      phi0 <- phi <- Phi <- 0
      if (d + D <= 1) {
        phi0 <- runif(1, -3, 3)
      }
      if (p[i] > 0) {
        phi <- stationary_ar(p[i])
      }
      if (seasonal_periods[i] > 1 & P > 0) {
        Phi <- stationary_ar(P)
      }
      pi_weights <- c(
        phi0,
        pi_coefficients(ar = phi, sar = Phi, d = d, D = D, m = seasonal_periods[i])
      )
      ar[seq_along(pi_weights), i] <- pi_weights
    }
  }
  # Trim redundant zeros
  nonzero <- which(rowSums(abs(ar)) > 0)
  ar <- ar[seq(max(nonzero)), , drop = FALSE]

  # Return mar object
  structure(list(ar = ar, sigmas = sigmas, weights = weights, m = seasonal_periods),
    class = "mar"
  )
}
