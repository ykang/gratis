#' @importFrom generics generate
#' @export
generics::generate

#' Generate a tsibble of synthetic data from a Mixture Autoregressive model
#'
#' This function simulates multiple random sample paths from a mixture of k Gaussian AR(p) processes.
#' The model is of the form
#' \deqn{y_t = \phi_{0,i} + \phi_{1,i}y_{t-1} + \dots + \phi_{p,i}y_{t-p} + \sigma_{i,t}\epsilon_t}
#' with probability \eqn{\alpha_i}, where \eqn{\epsilon_t} is a N(0,1) variate.
#' The index of the tsibble is guessed from the MAR model seasonal periods.
#' @param x A `mar` object, usually the output of \code{\link{mar_model}()}.
#' @param length length of series to generate
#' @param nseries number of series to generate
#' @param ... Other arguments, passed to \code{\link{simulate.mar}}.
#' @return `tsibble` object with `length` rows and 3 columns.
#' @references Feng Li, Mattias Villani, and Robert Kohn. (2010). Flexible Modeling of
#'     Conditional Distributions using Smooth Mixtures of Asymmetric Student T Densities,
#'     Journal of Statistical Planning and Inference, 140(12), pp. 3638-3654.
#' @author Rob J Hyndman
#' @seealso \code{\link{mar_model}}, \code{\link{simulate.mar}}
#' @examples
#' # MAR model with constant variances
#' ar <- cbind(c(0, 0.8, 0), c(0, 0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model1 <- mar_model(ar = ar, sigmas = c(1, 2), weights = weights)
#' generate(model1, nseries = 5)
#' # MAR model for hourly data with daily and weekly periods
#' hourly_model <- mar_model(seasonal_periods = c(24, 24*7))
#' generate(hourly_model)
#' @export
generate.mar <- function(x, length = 100, nseries = 10, ...) {
  generate_gratis(x, length, nseries, ...)
}

#' @rdname generate.mar
#' @export
generate.ets <- function(x, length = 100, nseries = 10, ...) {
  generate_gratis(x, length, nseries, ...)
}

#' @rdname generate.mar
#' @export
generate.Arima <- function(x, length = 100, nseries = 10, ...) {
  generate_gratis(x, length, nseries, ...)
}

# Generic generate function for models in gratis package
generate_gratis <- function(x, length = 100, nseries = 10, ...) {
  if (!is.null(x$m)) {
    m <- x$m
  } else if (!is.null(x$frequency)) {
    m <- x$frequency
  } else {
    m <- 1
  }
  tsmatrix <- forecast::msts(matrix(0, nrow = length, ncol = nseries), 
                seasonal.periods = unique(m), ts.frequency = max(m))
  for (i in seq(nseries)) {
    tsmatrix[, i] <- simulate(x, nsim = length, ...)
  }
  out <- tsibble::as_tsibble(tsmatrix)
  freq <- min(m)
  if(abs(freq - 365.25/7) < 1e-4 | freq == 52) {
    # Weekly data (which doesn't convert properly in tsibble)
    out$index <- seq(as.POSIXlt("01-01-01", tz="UTC"), by="7 days", length=NROW(out))
    out$index <- tsibble::yearweek(out$index)
  } else if(freq == 24) {
    # Hourly data
    out$index <- seq(as.POSIXlt("01-01-01", tz="UTC"), by="1 hour", length=NROW(out))
  } else if(freq == 48) {
    # Half hourly data
    out$index <- seq(as.POSIXlt("01-01-01", tz="UTC"), by="30 min", length=NROW(out))
  } else if(freq == 96) {
    # Quarter hour data
    out$index <- seq(as.POSIXlt("01-01-01", tz="UTC"), by="15 min", length=NROW(out))
  } else if(freq == 60) {
    # Minute data
    out$index <- seq(as.POSIXlt("01-01-01", tz="UTC"), by="1 min", length=NROW(out))
  }
  # Return
  tsibble::as_tsibble(out, index=index, key=key)
}
