#' Specify an ARIMA model
#'
#' Construct a Gaussian ARIMA\eqn{(p,d,q)(P,D,Q)_m} model that can be used with
#' \code{\link{simulate.gratis_arima}()} or \code{\link{generate.Arima}()}.
#' Omitted orders and coefficients are inferred where possible or randomly
#' selected from the supported parameter space. Construction does not simulate
#' observations or store a seed for later simulations.
#'
#' When orders are omitted, \code{p} and \code{q} are sampled from
#' \eqn{\{0,1,2,3\}}, \code{P} and \code{Q} from \eqn{\{0,1,2\}},
#' \code{d} from \eqn{\{0,1,2\}}, and \code{D} from \eqn{\{0,1\}}, with the
#' restriction \eqn{d + D \le 2}. If \code{constant} is omitted, it is sampled
#' from \eqn{(-3,3)} unless \eqn{d + D = 2}, in which case it is set to zero.
#' AR and MA parameters are sampled so the stationary and invertible parts of
#' the process are valid, and \code{sigma} is sampled from \eqn{(1,5)}.
#' Above order two, random AR polynomials allow every possible count of real
#' roots. This has full support, but is not uniform over coefficient space.
#'
#' @param frequency Positive integer seasonal period, for example 12 for monthly
#'   data. A value of 1 denotes a non-seasonal model.
#' @param p,q Non-negative integer non-seasonal AR and MA orders. If omitted,
#'   inferred from supplied \code{phi} and \code{theta}, respectively, or sampled.
#' @param d Non-negative integer ordinary difference order, with \code{d + D <= 2}.
#' @param P,Q Non-negative integer seasonal AR and MA orders. If omitted,
#'   inferred from supplied \code{Phi} and \code{Theta}, respectively, or sampled.
#' @param D Non-negative integer seasonal difference order. Seasonal arguments
#'   require \code{frequency > 1}.
#' @param constant Mean when \code{d + D = 0}, or slope of a linear drift
#'   when \code{d + D = 1}. Must be zero when \code{d + D > 1}.
#' @param phi Finite numeric vector of length \code{p} containing AR parameters.
#' @param theta Finite numeric vector of length \code{q} containing MA parameters.
#'   MA lag polynomials use positive signs: \eqn{1 + \theta(B)}.
#' @param Phi Finite numeric vector of length \code{P} containing seasonal AR
#'   parameters, at lags \code{frequency}, \code{2 * frequency}, and so on.
#' @param Theta Finite numeric vector of length \code{Q} containing seasonal MA
#'   parameters, using the same positive-sign convention as \code{theta}.
#' @param sigma The standard deviation of the noise. Its squared variance must
#'   lie in the positive normal floating-point range, in both R and Python.
#' @return A \code{gratis_arima} object inheriting from \code{Arima}. Its
#'   simulation method uses scaled innovations and preserves all seasonal
#'   coefficients, including coefficients whose sum is zero. Models fitted
#'   separately with \pkg{forecast} retain their own simulation method.
#' @author Rob J Hyndman
#' @seealso \code{\link{simulate.gratis_arima}} for initialization and controlled
#'   innovations; \code{\link{generate.Arima}} for batches of paths.
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
  if (length(frequency) != 1L || !is.finite(frequency) ||
      frequency < 1 || frequency != as.integer(frequency)) {
    stop("frequency must be a positive integer")
  }
  # sigma
  if (is.null(sigma)) {
    # Choose variance in [1,5]
    sigma <- runif(1, 1, 5)
  }
  validate_sigma(sigma)
  # Model orders
  if (is.null(p)) {
    p <- if (is.null(phi)) sample(c(0, 1, 2, 3), 1) else length(phi)
  }
  if (is.null(d)) {
    d <- sample(c(0, 1, 2), 1)
  }
  if (is.null(q)) {
    q <- if (is.null(theta)) sample(c(0, 1, 2, 3), 1) else length(theta)
  }
  if (frequency > 1) {
    if (is.null(P)) {
      P <- if (is.null(Phi)) sample(c(0, 1, 2), 1) else length(Phi)
    }
    if (is.null(D)) {
      if (d == 2) {
        D <- 0
      } else {
        D <- sample(c(0, 1), 1)
      }
    }
    if (is.null(Q)) {
      Q <- if (is.null(Theta)) sample(c(0, 1, 2), 1) else length(Theta)
    }
  } else {
    if (any(c(P, D, Q) != 0) || length(Phi) > 0L || length(Theta) > 0L) {
      stop("seasonal arguments require frequency greater than 1")
    }
    P <- D <- Q <- 0
  }
  orders <- list(p = p, d = d, q = q, P = P, D = D, Q = Q)
  for (name in names(orders)) {
    value <- orders[[name]]
    if (length(value) != 1L || !is.finite(value) || value < 0 ||
        value != as.integer(value)) {
      stop(paste(name, "must be a non-negative integer"))
    }
  }
  if (d + D > 2) stop("d + D must be no greater than 2")
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
  if (length(c) != 1L || !is.finite(c)) stop("constant must be finite")
  if (d + D > 1 && c != 0) stop("constant must be zero when d + D > 1")

  # Now put it in an Arima class by "fitting" the model to a ts object
  pars <- c(phi, theta, Phi, Theta)
  if (c != 0) {
    pars <- c(pars, c)
  }
  n <- max(100, (d + D + 1) * frequency)
  model <- forecast::Arima(ts(rep(0, n), frequency = frequency, end = 0),
    order = c(p, d, q), seasonal = c(P, D, Q), fixed = pars, include.constant = (c != 0)
  )
  model$x <- model$fitted <- model$residuals <- model$series <- model$call <- NULL
  model$nobs <- model$n.cond <- model$aicc <- model$bic <- model$var.coef <- NULL
  model$mask <- model$code <- NULL
  model$sigma2 <- sigma^2
  model$frequency <- frequency
  model$loglik <- model$aic <- NA
  class(model) <- c("gratis_arima", "forecast_ARIMA", "Arima", "ARIMA")

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
    counts <- seq(p %% 2, p, by = 2)
    n_real_roots <- counts[sample.int(length(counts), 1L)]
    inv_real_roots <- runif(n_real_roots, -1, 1)
    # Generate inverse complex roots in conjugate pairs
    n_complex_roots <- (p - n_real_roots) / 2
    r <- runif(n_complex_roots, -1, 1)
    angle <- runif(n_complex_roots, -pi, pi)
    inv_complex_roots <- c(complex(argument = angle, modulus = r),
                           complex(argument = -angle, modulus = r))
    # Multiply (1 - inverse_root * B), without dividing by near-zero roots.
    poly <- 1
    for (root in c(inv_real_roots, inv_complex_roots)) {
      poly <- c(poly, 0) - c(0, root * poly)
    }
    phi <- -Re(poly[-1L])
  }
  return(phi)
}

validate_sigma <- function(sigma) {
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be positive and finite")
  }
  if (sigma < sqrt(.Machine$double.xmin) || sigma > sqrt(.Machine$double.xmax)) {
    stop("sigma squared must be a positive finite normal variance")
  }
  invisible(NULL)
}

#' Compute AR recursion coefficients from SARIMA coefficients
#'
#' Convert SARIMA coefficients into the equivalent infinite-order AR recursion
#' coefficients, truncated after the last coefficient whose absolute value is
#' greater than \code{tol}. These coefficients are used internally by
#' \code{\link{mar_model}()} and are also exported for users who need the same
#' conversion.
#'
#' AR lag polynomials use negative coefficient signs, while MA lag polynomials
#' use positive signs. The result expands the ratio of the AR/differencing
#' polynomial to the MA polynomial; no stationarity or invertibility is imposed.
#'
#' @param ar,ma Finite numeric vectors of non-seasonal AR and MA coefficients,
#'   excluding the leading 1 of each lag polynomial.
#' @param sar,sma Finite numeric vectors of seasonal AR and MA coefficients.
#'   Ignored when \code{m = 1}.
#' @param d,D Non-negative integer ordinary and seasonal difference orders.
#'   Seasonal differencing is ignored when \code{m = 1}.
#' @param m Positive integer seasonal period. A value of 1 is non-seasonal.
#' @param tol Non-negative finite truncation tolerance. Coefficients are returned through
#'   the last position whose absolute value exceeds \code{tol}.
#' @param n_terms Lag limit for an infinite AR approximation with MA terms.
#'   The limit is increased to include the full AR and MA polynomial degrees.
#'   A warning is issued if the tail has not decayed below \code{tol}.
#'   Without MA terms the complete finite AR polynomial is returned.
#'
#' @return A numeric vector of AR recursion coefficients in increasing lag
#'   order, without an intercept or leading 1. Complete cancellation returns
#'   \code{numeric(0)}. An unconverged MA approximation issues a warning;
#'   increase \code{n_terms}, noting that noninvertible MA terms may not converge.
#' @seealso \code{\link{mar_model}}
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
pi_coefficients <- function(ar = 0, d = 0L, ma = 0, sar = 0, D = 0L, sma = 0, m = 1L, tol = 1e-07, n_terms = 500L) {
  if (any(!is.finite(c(ar, ma, sar, sma)))) stop("AR and MA coefficients must be finite")
  for (value in list(d, D, m)) {
    if (length(value) != 1L || !is.finite(value) || value < 0 || value != floor(value)) {
      stop("d, D and m must be non-negative integers")
    }
  }
  if (m < 1L) stop("m must be positive")
  if (length(n_terms) != 1L || !is.finite(n_terms) || n_terms < 1 ||
      n_terms != as.integer(n_terms)) stop("n_terms must be a positive integer")
  if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("tol must be non-negative and finite")
  }
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

  numerator <- as.numeric(stats::coef(ar * sar))
  denominator <- as.numeric(stats::coef(ma * sma))
  q <- length(denominator) - 1L
  if (q == 0L || all(denominator[-1L] == 0)) {
    coefficients <- -numerator[-1L]
  } else {
    # Series division A(B) / M(B), using the standard positive MA convention.
    n <- max(n_terms, length(numerator) - 1L, q)
    pie <- numeric(n + 1L)
    pie[1L] <- 1
    for (j in seq_len(n)) {
      lags <- seq_len(min(j, q))
      value <- if (j < length(numerator)) numerator[j + 1L] else 0
      pie[j + 1L] <- value - sum(denominator[lags + 1L] * pie[j - lags + 1L])
    }
    coefficients <- -pie[-1L]
    if (any(!is.finite(coefficients)) || any(abs(utils::tail(coefficients, q)) > tol)) {
      warning("MA approximation has not decayed below tol; increase n_terms")
    }
  }
  nonzero <- which(abs(coefficients) > tol)
  if (length(nonzero) == 0L) return(numeric(0))
  head(coefficients, max(nonzero))
}
