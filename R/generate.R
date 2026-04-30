#' @importFrom generics generate
#' @export
generics::generate

#' Generate a tsibble of synthetic time series
#'
#' Generate one or more synthetic time series from a fitted or specified model.
#' Methods are provided for \code{\link{mar_model}()} objects, ETS objects from
#' \code{\link{ets_model}()}, and ARIMA objects from \code{\link{arima_model}()}.
#' Each series is produced by repeatedly calling the corresponding
#' \code{\link[stats]{simulate}()} method.
#'
#' The generated series are returned as a tsibble. The index is inferred from
#' the model frequency or seasonal periods where possible, with special handling
#' for common weekly, hourly, half-hourly, quarter-hourly, and minute data.
#'
#' @param x A model object of class \code{mar}, \code{ets}, or \code{Arima}.
#' @param length Length of each series to generate.
#' @param nseries Number of series to generate.
#' @param ... Other arguments passed to the relevant \code{simulate()} method.
#'   For \code{mar} objects these are passed to \code{\link{simulate.mar}()}.
#' @return A tsibble containing the generated values. For multiple series the
#'   output has \code{length * nseries} rows and includes a key column identifying
#'   each series.
#' @references Feng Li, Mattias Villani, and Robert Kohn. (2010). Flexible Modeling of
#'     Conditional Distributions using Smooth Mixtures of Asymmetric Student T Densities,
#'     Journal of Statistical Planning and Inference, 140(12), pp. 3638-3654.
#' @author Rob J Hyndman
#' @seealso \code{\link{mar_model}}, \code{\link{arima_model}},
#'   \code{\link{ets_model}}, \code{\link{simulate.mar}}
#' @examples
#' # MAR model with constant variances
#' phi <- cbind(c(0.8, 0), c(0.6, 0.3))
#' weights <- c(0.8, 0.2)
#' model1 <- mar_model(phi = phi, sigmas = c(1, 2), weights = weights)
#' generate(model1, nseries = 5)
#'
#' # MAR model for hourly data with daily and weekly periods
#' hourly_model <- mar_model(seasonal_periods = c(24, 24*7))
#' generate(hourly_model, length = 24 * 7, nseries = 2)
#'
#' # ARIMA and ETS models also work
#' generate(arima_model(p = 1, d = 0, q = 0, phi = 0.5), length = 20, nseries = 2)
#' generate(ets_model(error = "A", trend = "N", seasonal = "N"), length = 20, nseries = 2)
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
  if("ARIMA" %in% class(x)) {
    m <- x$arma[5]
  } else if (!is.null(x$m)) {
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
  make_tsibble(tsmatrix, m)
}

# Create tsibble from ts matrix
make_tsibble <- function(tsmatrix, seasonal_periods) {
  if(!("ts" %in% class(tsmatrix)))
    tsmatrix <- ts(tsmatrix)
  out <- tsibble::as_tsibble(tsmatrix)
  if("key" %in% colnames(out))
    nseries <- length(unique(out$key))
  else
    nseries <- 1
  length <- NROW(out)/nseries
  freq <- min(seasonal_periods)
  if(abs(freq - 365.25/7) < 1e-4 | freq == 52) {
    # Weekly data (which doesn't convert properly in tsibble)
    out$index <- rep(seq(as.POSIXlt("01-01-01", tz="UTC"), by="7 days", length=length), nseries)
    out$index <- tsibble::yearweek(out$index)
  } else if(freq == 24) {
    # Hourly data
    out$index <- rep(seq(as.POSIXlt("01-01-01", tz="UTC"), by="1 hour", length=length), nseries)
  } else if(freq == 48) {
    # Half hourly data
    out$index <- rep(seq(as.POSIXlt("01-01-01", tz="UTC"), by="30 min", length=length), nseries)
  } else if(freq == 96) {
    # Quarter hour data
    out$index <- rep(seq(as.POSIXlt("01-01-01", tz="UTC"), by="15 min", length=length), nseries)
  } else if(freq == 60) {
    # Minute data
    out$index <- rep(seq(as.POSIXlt("01-01-01", tz="UTC"), by="1 min", length=length), nseries)
  }
  # Return
  if("key" %in% colnames(out))
    return(tsibble::as_tsibble(out, index=index, key=key))
  else
    return(tsibble::as_tsibble(out, index=index))
}

utils::globalVariables(c("index","key"))
