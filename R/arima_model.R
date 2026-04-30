#' Specify an ARIMA model
#'
#' Construct a Gaussian ARIMA\eqn{(p,d,q)(P,D,Q)_m} model that can be used with
#' \code{\link[forecast]{simulate.Arima}()} or \code{\link{generate.Arima}()}.
#' Any omitted argument is randomly selected from the supported parameter space.
#'
#' When orders are omitted, \code{p} and \code{q} are sampled from
#' \eqn{\{0,1,2,3\}}, \code{P} and \code{Q} from \eqn{\{0,1,2\}},
#' \code{d} from \eqn{\{0,1,2\}}, and \code{D} from \eqn{\{0,1\}}, with the
#' restriction \eqn{d + D \le 2}. If \code{constant} is omitted, it is sampled
#' from \eqn{(-3,3)} unless \eqn{d + D = 2}, in which case it is set to zero.
#' AR and MA parameters are sampled so the stationary and invertible parts of
#' the process are valid, and \code{sigma} is sampled from \eqn{(1,5)}.
#'
#' @param frequency The length of the seasonal period (e.g., 12 for monthly data).
#' @param p Non-seasonal autoregressive order.
#' @param d Non-seasonal differencing order.
#' @param q Non-seasonal moving-average order.
#' @param P Seasonal autoregressive order.
#' @param D Seasonal differencing order.
#' @param Q Seasonal moving-average order.
#' @param constant Intercept term.
#' @param phi A numeric p-vector containing the AR parameters.
#' @param theta A numeric q-vector containing the MA parameters.
#' @param Phi A numeric P-vector containing the seasonal AR parameters.
#' @param Theta A numeric Q-vector containing the seasonal MA parameters.
#' @param sigma The standard deviation of the noise.
#' @return An \code{Arima} object compatible with the \pkg{forecast} package's
#'   simulation methods.
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
#' plot(simulate(model1, 100))
#' plot(simulate(model2, 100))
#' plot(simulate(model3, 100))
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

#' Compute AR recursion coefficients from SARIMA coefficients
#'
#' Convert SARIMA coefficients into the equivalent infinite-order AR recursion
#' coefficients, truncated after the last coefficient whose absolute value is
#' greater than \code{tol}. These coefficients are used internally by
#' \code{\link{mar_model}()} and are also exported for users who need the same
#' conversion.
#'
#' @param ar Non-seasonal AR coefficients.
#' @param d Number of non-seasonal differences.
#' @param ma Non-seasonal MA coefficients.
#' @param sar Seasonal AR coefficients.
#' @param D Number of seasonal differences.
#' @param sma Seasonal MA coefficients.
#' @param m Seasonal period.
#' @param tol Tolerance used for truncation. Coefficients are returned through
#'   the last position whose absolute value exceeds \code{tol}.
#'
#' @return A numeric vector of AR recursion coefficients.
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
#' pi_coefficients(ar = 0.8)
#' pi_coefficients(ar = c(0.4, -0.2), d = 1)
#' pi_coefficients(ar = 0.3, sar = 0.2, D = 1, m = 12)
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
  theta <- -c(stats::coef(ma * sma))[-1]
  if (length(theta) == 0L) {
    theta <- 0
  }
  phi <- -c(stats::coef(ar * sar)[-1], numeric(n))
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
