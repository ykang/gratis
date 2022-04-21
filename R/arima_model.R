#' Specify parameters for an ARIMA model
#'
#' This function allows the parameters of a Gaussian \eqn{ARIMA(p,d,q)(P,D,Q)[m]}
#' process to be specified. The output can be used in \code{\link[forecast]{simulate.Arima}()}
#' and \code{\link{generate.Arima}}.
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
#' @param p An integer equal to the non-seasonal autoregressive order
#' @param d An integer equal to the non-seasonal order of differencing
#' @param q An integer equal to the non-seasonal moving average order
#' @param P An integer equal to the seasonal autoregressive order
#' @param D An integer equal to the seasonal order of differencing
#' @param Q An integer equal to the seasonal moving average order
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
#' model1 <- arima_model(p = 2, d = 0, q = 0)
#' # An AR(2) model with specific parameters
#' model2 <- arima_model(p = 2, d = 0, q = 0, phi = c(1.34, -0.64), sigma = 15)
#' # Seasonal ARIMA model with randomly selected parameters
#' model3 <- arima_model(frequency = 4)
#' # Simulate from each model and plot the results
#' library(forecast)
#' simulate(model1, 100) %>% plot()
#' simulate(model2, 100) %>% plot()
#' simulate(model3, 100) %>% plot()
#' @export
arima_model <- function(frequency = 1, p = NULL, d = NULL, q = NULL,
                        P = NULL, D = NULL, Q = NULL, constant = NULL,
                        phi = NULL, theta = NULL, Phi = NULL, Theta = NULL,
                        sigma = NULL) {
  # sigma
  if (is.null(sigma)) {
    # Choose variance in [1,5]
    sigma <- runif(1, 1, 5)
  }
  # Model orders
  if (is.null(p)) {
    p <- sample(c(0, 1, 2, 3), 1)
  }
  if (is.null(d)) {
    d <- sample(c(0, 1, 2), 1)
  }
  if (is.null(q)) {
    q <- sample(c(0, 1, 2, 3), 1)
  }
  if (frequency > 1) {
    if (is.null(P)) {
      P <- sample(c(0, 1, 2), 1)
    }
    if (is.null(D)) {
      if (d == 2) {
        D <- 0
      } else {
        D <- sample(c(0, 1), 1)
      }
    }
    if (is.null(Q)) {
      Q <- sample(c(0, 1, 2), 1)
    }
  } else {
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
  class(model) <- c("forecast_ARIMA", "Arima", "ARIMA")

  return(model)
}

# Function to return random coefficients from a stationary AR(p) process
stationary_ar <- function(p) {
  p <- as.integer(p)
  if (p < 1) {
    stop("p must be a positive integer")
  } else if (p == 1L) {
    phi <- runif(1, -1, 1)
  } else if (p == 2L) {
    phi2 <- runif(1, -1, 1)
    phi1 <- runif(1, phi2 - 1, 1 - phi2)
    phi <- c(phi1, phi2)
  } else {
    # Generate inverse real roots
    n_real_roots <- p %% 2
    inv_real_roots <- runif(n_real_roots, -1, 1)
    # Generate inverse complex roots in conjugate pairs
    n_complex_roots <- (p - n_real_roots) / 2
    r <- runif(n_complex_roots, -1, 1)
    angle <- runif(n_complex_roots, -pi, pi)
    inv_complex_roots <- c(complex(argument = angle, modulus = r),
                           complex(argument = -angle, modulus = r))
    # Find polynomial with these as roots
    poly <- suppressWarnings(polynom::poly.calc(1/c(inv_real_roots, inv_complex_roots)))
    # Scale to have constant 1
    poly <- poly / poly[1]
    phi <- -as.numeric(poly)[-1]
  }
  return(phi)
}

#' Compute pi coefficients of an AR process from SARIMA coefficients.
#'
#' Convert SARIMA coefficients to pi coefficients of an AR process.
#' @param ar AR coefficients in the SARIMA model.
#' @param d number of differences in the SARIMA model.
#' @param ma MA coefficients in the SARIMA model.
#' @param sar seasonal AR coefficients in the SARIMA model.
#' @param D number of seasonal differences in the SARIMA model.
#' @param sma seasonal MA coefficients in the SARIMA model.
#' @param m seasonal period in the SARIMA model.
#' @param tol tolerance value used. Only return up to last element greater than tolerance.
#'
#' @return A vector of AR coefficients.
#' @author Rob J Hyndman
#' @export
#' @importFrom stats dist
#' @importFrom stats tsp
#' @importFrom stats tsp<-
#' @importFrom stats acf
#' @importFrom forecast BoxCox.lambda
#' @importFrom stats loess
#'
#' @examples
#' # Not Run
pi_coefficients <- function(ar = 0, d = 0L, ma = 0, sar = 0, D = 0L, sma = 0, m = 1L, tol = 1e-07) {
  # non-seasonal AR
  ar <- polynomial(c(1, -ar)) * polynomial(c(1, -1))^d

  # seasonal AR
  if (m > 1) {
    P <- length(sar)
    seasonal_poly <- numeric(m * P)
    seasonal_poly[m * seq(P)] <- sar
    sar <- polynomial(c(1, -seasonal_poly)) * polynomial(c(1, rep(0, m - 1), -1))^D
  } else {
    sar <- 1
  }

  # non-seasonal MA
  ma <- polynomial(c(1, ma))

  # seasonal MA
  if (m > 1) {
    Q <- length(sma)
    seasonal_poly <- numeric(m * Q)
    seasonal_poly[m * seq(Q)] <- sma
    sma <- polynomial(c(1, seasonal_poly))
  } else {
    sma <- 1
  }

  n <- 500L
  theta <- -c(coef(ma * sma))[-1]
  if (length(theta) == 0L) {
    theta <- 0
  }
  phi <- -c(coef(ar * sar)[-1], numeric(n))
  q <- length(theta)
  pie <- c(numeric(q), 1, numeric(n))
  for (j in seq(n)) {
    pie[j + q + 1L] <- -phi[j] - sum(theta * pie[(q:1L) + j])
  }
  pie <- pie[(0L:n) + q + 1L]

  # Return up to last element greater than tol
  maxj <- max(which(abs(pie) > tol))
  pie <- head(pie, maxj)
  return(-pie[-1])
}
