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
  }
  else {
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
  }
  else {
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
  for (j in seq(n))
    pie[j + q + 1L] <- -phi[j] - sum(theta * pie[(q:1L) + j])
  pie <- pie[(0L:n) + q + 1L]

  # Return up to last element greater than tol
  maxj <- max(which(abs(pie) > tol))
  pie <- head(pie, maxj)
  return(-pie[-1])
}

#' Compute pi coefficients from ARIMA model
#'
#' Compute pi coefficients from ARIMA model
#' @param object An object of class "Arima"
#'
#' @return A vector of AR coefficients
#' @author Rob J Hyndman
#' @export
#'
#' @examples
#' # Not Run
arinf <- function(object) {
  if (!("Arima" %in% class(object))) {
    stop("Argument should be an ARIMA object")
  }
  pi_coefficients(
    ar = object$model$phi, ma = object$model$theta,
    d = object$arma[6], D = object$arma[7], m = object$arma[5]
  )
}

# library(forecast)
# USAccDeaths %>% auto.arima %>% arinf %>% plot
# lynx %>% auto.arima %>% arinf %>% plot

#' Set the number of seasonal differences for yearly data to be -1.
#'
#' @param x Univariate time series or numerical vector
#'
#' @return
#' NA
#'
#' @examples
#' # Not Run
nsdiffs1 <- function(x) {
  c(nsdiffs = ifelse(frequency(x) == 1L, -1, forecast::nsdiffs(x)))
}
