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
#' Automatic burn-in for models created by \code{\link{arima_model}()} is
#' calculated once per batch. This does not change the simulation behavior of
#' separately fitted \pkg{forecast} models.
#'
#' The generated series are returned as a tsibble. The index is inferred from
#' the model frequency or seasonal periods where possible, with special handling
#' for common weekly, hourly, half-hourly, quarter-hourly, and minute data.
#'
#' @param x A model object of class \code{mar}, \code{ets}, or \code{Arima}.
#' @param length Positive integer output length of each series, excluding burn-in.
#' @param nseries Positive integer number of series to generate.
#' @param seed Optional seed set once for the whole batch. Successive series
#'   use successive random draws, and the previous RNG state is restored after
#'   generation. The default, \code{NULL}, advances the current RNG normally.
#' @param ... Other arguments passed to the relevant \code{simulate()} method.
#'   See \code{\link{simulate.mar}()} and \code{\link{simulate.gratis_arima}()}
#'   for initialization controls. Supplied innovations must produce exactly
#'   \code{length} observations. Not every model accepts the same arguments.
#' @return A tsibble containing the generated values. For multiple series the
#'   output has \code{length * nseries} rows and includes a key column identifying
#'   each series. Unlike the Python API, the result contains calendar or seasonal
#'   index information rather than being an unlabelled numeric matrix.
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
generate.mar <- function(x, length = 100, nseries = 10, ..., seed = NULL) {
  generate_gratis(x, length, nseries, ..., seed = seed)
}

#' @rdname generate.mar
#' @export
generate.ets <- function(x, length = 100, nseries = 10, ..., seed = NULL) {
  generate_gratis(x, length, nseries, ..., seed = seed)
}

#' @rdname generate.mar
#' @export
generate.Arima <- function(x, length = 100, nseries = 10, ..., seed = NULL) {
  generate_gratis(x, length, nseries, ..., seed = seed)
}

# Generic generate function for models in gratis package
generate_gratis <- function(x, length = 100, nseries = 10, ..., seed = NULL) {
  for (value in list(length, nseries)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      stop("length and nseries must be positive integers")
    }
  }
  simulation_args <- list(...)
  if (inherits(x, "gratis_arima") && is.null(simulation_args$n.start) &&
      is.null(simulation_args$start.innov)) {
    simulation_args$n.start <- arima_burnin(x)
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
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
    tsmatrix[, i] <- do.call(simulate, c(list(object = x, nsim = length), simulation_args))
  }
  make_tsibble(tsmatrix, m)
}

# Create tsibble from ts matrix
make_tsibble <- function(tsmatrix, seasonal_periods) {
  if(!("ts" %in% class(tsmatrix)))
    tsmatrix <- forecast::msts(tsmatrix, seasonal.periods = unique(seasonal_periods),
                             ts.frequency = max(seasonal_periods))
  freq <- min(seasonal_periods)
  weekly <- abs(freq - 365.25/7) < 1e-4 || freq == 52
  # Weekly indices are assigned explicitly below; avoid a conflicting tsibble
  # conversion of a 52-observation calendar year.
  out <- tsibble::as_tsibble(if (weekly) stats::ts(tsmatrix, frequency = 1) else tsmatrix)
  if("key" %in% colnames(out))
    nseries <- length(unique(out$key))
  else
    nseries <- 1
  length <- NROW(out)/nseries
  if(weekly) {
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
