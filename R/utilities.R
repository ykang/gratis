
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
    lambda <- forecast::BoxCox.lambda(x, lower = 0, upper = 1)
    bc.x <- as.numeric(forecast::BoxCox(x, lambda))
    seasonal <- stats::stl(ts(bc.x, frequency = ppy), s.window = "periodic")$time.series[, 1]
    bc.x <- bc.x - as.numeric(seasonal)
    x <- as.numeric(forecast::InvBoxCox(bc.x, lambda)) + x - x
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
