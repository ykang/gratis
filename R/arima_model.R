#' Specify parameters for an ARIMA model
#'
#' This function allows the parameters of a Gaussian \eqn{ARIMA(p,d,q)(P,D,Q)[m]}
#' process to be specified. The output can be used in \code{\link[forecast]{simulate.Arima}()}.
#' If any argument is \code{NULL}, the corresponding parameters are randomly selected.
#' The AR and MA orders p and q are chosen from \{0,1,2,3\}, the seasonal AR and MA
#' orders P and Q are from \{0,1,2\}, while the order of differencing,
#' d is in \{0,1,2\}, and the order of seasonal differencing D is in \{0,1\}, with the
#' restriction that \eqn{d+D \le 2}. If \code{constant} is \code{NULL}, it is set to
#' 0 if \eqn{d+D = 2}, otherwise it is uniformly sampled on (-3,3).
#' The model orders and the parameters are uniformly sampled. The AR and MA parameters are selected
#' to give stationary and invertible processes when \eqn{d=D=0}. The noise variance sigma
#' is uniformly sampled on (1,5). The parameterization is as specified in Hyndman & Athanasopoulos (2021).
#'
#' @param frequency The length of the seasonal period (e.g., 12 for monthly data).
#' @param order A numeric vector specifying the non-seasonal part of the ARIMA model: the three components \code{(p, d, q)}.
#' @param seasonal A numeric vector specifying the seasonal part of the ARIMA model: the three components \code{(P, D, Q)}.
#' @param constant The intercept term
#' @param phi A numeric p-vector containing the AR parameters.
#' @param theta A numeric p-vector containing the MA parameters.
#' @param Phi A numeric p-vector containing the seasonal AR parameters.
#' @param Theta A numeric p-vector containing the seasonal MA parameters.
#' @param sigma The standard deviation of the noise.
#' @return An `Arima` object as described in the \code{\link[stats]{arima}} function from the stats package.
#' @author Rob J Hyndman
#' @seealso \code{\link[forecast]{simulate.Arima}}
#' @examples
#' # An AR(2) model with random parameters
#' model1 <- arima_model(order = c(2, 0, 0))
#' # An AR(2) model with specific parameters
#' model2 <- arima_model(order = c(2, 0, 0), phi = c(1.34, -0.64), sigma = 15)
#' # Seasonal ARIMA model with randomly selected parameters
#' model3 <- arima_model(frequency = 4)
#' # Simulate from each model and plot the results
#' library(forecast)
#' simulate(model1, 100) %>% plot()
#' simulate(model2, 100) %>% plot()
#' simulate(model3, 100) %>% plot()
#' @export
arima_model <- function(frequency = 1, order = NULL, seasonal = NULL, constant = NULL,
                        phi = NULL, theta = NULL, Phi = NULL, Theta = NULL,
                        sigma = NULL) {
  # sigma
  if (is.null(sigma)) {
    # Choose variance in [1,5]
    sigma <- runif(1, 1, 5)
  }
  # Model orders
  if (!is.null(order)) {
    p <- order[1]
    d <- order[2]
    q <- order[3]
  } else {
    p <- sample(c(0, 1, 2, 3), 1)
    d <- sample(c(0, 1, 2), 1)
    q <- sample(c(0, 1, 2, 3), 1)
  }
  if (frequency > 1) {
    if (!is.null(seasonal)) {
      P <- seasonal[1]
      D <- seasonal[2]
      Q <- seasonal[3]
    } else {
      P <- sample(c(0, 1, 2), 1)
      if (d == 2) {
        D <- 0
      } else {
        D <- sample(c(0, 1), 1)
      }
      Q <- sample(c(0, 1, 2), 1)
    }
  } else {
    if (!is.null(seasonal)) {
      stop("frequency must be greater than 1 for seasonal models")
    }
    P <- D <- Q <- 0
  }
  # phi
  if (!is.null(phi)) {
    if (length(phi) != p) {
      stop("Dimension of phi does not match model order")
    }
  } else if (p > 0) {
    phi <- stationary_ar(p)
  }
  # theta
  if (!is.null(theta)) {
    if (length(theta) != q) {
      stop("Dimension of theta does not match model order")
    }
  } else if (q > 0) {
    theta <- -stationary_ar(q)
  }
  # Phi
  if (!is.null(Phi)) {
    if (length(Phi) != P) {
      stop("Dimension of Phi does not match model order")
    }
  } else if (P > 0) {
    Phi <- stationary_ar(P)
  }
  # Theta
  if (!is.null(Theta)) {
    if (length(Theta) != Q) {
      stop("Dimension of Theta does not match model order")
    }
  } else if (Q > 0) {
    Theta <- -stationary_ar(Q)
  }
  # Intercept
  if (!is.null(constant)) {
    c <- constant
  } else if (d + D <= 1) {
    c <- runif(1, -3, 3)
  } else {
    c <- 0
  }

  # Now put it in an Arima class by "fitting" the model to a ts object
  pars <- c(phi, theta, Phi, Theta)
  if (c != 0) {
    pars <- c(pars, c)
  }
  model <- forecast::Arima(ts(rep(0, 100), frequency = frequency, end = 0),
    order = c(p, d, q), seasonal = c(P, D, Q), fixed = pars, include.constant = (c != 0)
  )
  model$x <- model$fitted <- model$residuals <- model$series <- model$call <- NULL
  model$nobs <- model$n.cond <- model$aicc <- model$bic <- model$var.coef <- NULL
  model$mask <- model$code <- NULL
  model$sigma2 <- sigma^2
  model$frequency <- frequency
  model$loglik <- model$aic <- NA

  return(model)
}
