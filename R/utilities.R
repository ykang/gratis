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

SeasonalityTest <- function(input, ppy) {
  if (length(input) < 3 * ppy) {
    test_seasonal <- FALSE
  } else {
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- 1.645 / sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- (abs(xacf[ppy]) > clim[ppy])
    if (is.na(test_seasonal)) {
      test_seasonal <- FALSE
    }
  }
  return(test_seasonal)
}

Smoothing_ts2 <- function(x, spanw, fh) {
  ppy <- frequency(x)
  trend <- seq_len(x)
  if (ppy > 1) {
    ST <- SeasonalityTest(x, ppy)
  } else {
    ST <- FALSE
  }
  if (ST) {
    lambda <- BoxCox.lambda(x, lower = 0, upper = 1)
    bc.x <- as.numeric(BoxCox(x, lambda))
    seasonal <- stl(ts(bc.x, frequency = ppy), s.window = "periodic")$time.series[, 1]
    bc.x <- bc.x - as.numeric(seasonal)
    x <- as.numeric(InvBoxCox(bc.x, lambda)) + x - x
    suppressWarnings(x.loess <- loess(x ~ trend, span = spanw / length(x), degree = 1))
    x <- as.numeric(x.loess$fitted) + x - x
    SIin <- seasonal
    SIout <- head(rep(seasonal[(length(seasonal) - ppy + 1):length(seasonal)], fh), fh)
  } else {
    suppressWarnings(x.loess <- loess(x ~ trend, span = spanw / length(x), degree = 1))
    x <- as.numeric(x.loess$fitted) + x - x
    SIin <- rep(0, length(x))
    SIout <- rep(0, fh)
    lambda <- 1
  }
  output <- list(series = x, seasonalIn = SIin, seasonal = SIout, lambda = lambda)
  return(output)
}
