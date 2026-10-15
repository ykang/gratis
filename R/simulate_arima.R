#' Simulate from a specified ARIMA model
#'
#' Simulate the stationary ARMA part from zero states, discard burn-in, and
#' then apply differencing inverses and the model's mean or linear drift.
#' All random innovations, including burn-in, have standard deviation
#' \code{sqrt(object$sigma2)}. This method applies only to models constructed
#' by \code{\link{arima_model}()}, not to fitted \pkg{forecast} models.
#'
#' @param object A model created by \code{\link{arima_model}()}.
#' @param nsim Non-negative integer output length, excluding burn-in. Supplied
#'   \code{innov} overrides the length with \code{length(innov)}.
#' @param seed Optional integer random seed. A supplied seed temporarily sets
#'   the RNG state, restoring the previous state or its absence afterwards.
#'   \code{NULL} uses and advances the current stream for unsupplied innovations.
#' @param n.start Non-negative integer burn-in length. \code{NULL} selects an
#'   automatic length; zero starts from zero states without burn-in. An explicit
#'   value skips the energy search, but not the AR stationarity check.
#' @param innov Optional finite numeric vector of output-period innovations,
#'   already on the model's noise scale. Values are not multiplied by
#'   \code{sqrt(object$sigma2)} again.
#' @param start.innov Optional finite numeric vector of burn-in innovations.
#'   Its length determines \code{n.start} when that argument is \code{NULL}.
#' @param ... Other arguments, not supported.
#' @return A \code{\link[stats]{ts}} object with the model's frequency and the
#'   requested output length, or \code{numeric(0)} for zero output length.
#' @section Automatic initialization:
#' For an AR model, impulse-response windows double from a power of two covering
#' at least 64 steps, \code{2 * (AR order + MA order + 1)}, and twelve time
#' constants of the slowest pole. Burn-in is half the first window whose latter
#' half contains at most \code{1e-8} of its energy. Repeated and seasonal poles
#' and near AR-MA cancellation are included. A pure MA model uses its expanded
#' MA order as burn-in. The same rule is used in Python.
#'
#' This is a numerical tail estimate, not exact stationary initialization.
#' Searches exceeding \code{2^22} steps stop with an error requesting an explicit
#' \code{n.start}. Differencing inverses and the mean or drift are applied only
#' after the stationary ARMA burn-in has been discarded.
#' @seealso \code{\link{arima_model}}, \code{\link{generate.Arima}}
#' @export
#' @examples
#' model <- arima_model(p = 1, d = 0, q = 0, phi = 0.99, constant = 0, sigma = 0.01)
#' simulate(model, nsim = 20, seed = 1)
#' # Supplied innovations are already scaled; zero burn-in gives zero prehistory.
#' impulse_model <- arima_model(p = 1, d = 0, q = 0, phi = 0.5, constant = 0, sigma = 1)
#' as.numeric(simulate(impulse_model, innov = c(1, 0, 0), n.start = 0))
simulate.gratis_arima <- function(object, nsim = 100, seed = NULL,
                                  n.start = NULL, innov = NULL,
                                  start.innov = NULL, ...) {
  if (length(list(...))) stop("unsupported ARIMA simulation arguments")
  check_innov <- function(x, name) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop(paste(name, "must be a finite numeric vector"))
    }
    as.numeric(x)
  }
  if (!is.null(innov)) {
    innov <- check_innov(innov, "innov")
    nsim <- length(innov)
  }
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 0 || nsim != as.integer(nsim)) {
    stop("nsim must be a non-negative integer")
  }
  if (nsim == 0L) return(numeric(0))
  ar <- object$model$phi
  ma <- object$model$theta
  ar <- head(ar, max(c(0L, which(ar != 0))))
  ma <- head(ma, max(c(0L, which(ma != 0))))
  if (!is.null(start.innov)) {
    start.innov <- check_innov(start.innov, "start.innov")
    if (is.null(n.start)) n.start <- length(start.innov)
  }
  if (is.null(n.start)) n.start <- arima_burnin(object) else arima_check_stationarity(object)
  if (length(n.start) != 1L || !is.finite(n.start) || n.start < 0 ||
      n.start != as.integer(n.start)) stop("n.start must be a non-negative integer")
  if (!is.null(start.innov) && length(start.innov) != n.start) {
    stop("start.innov must have length n.start")
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(seed)
  }
  sigma <- sqrt(object$sigma2)
  validate_sigma(sigma)
  if (is.null(start.innov)) start.innov <- rnorm(n.start, sd = sigma)
  if (is.null(innov)) innov <- rnorm(nsim, sd = sigma)

  # Zero padding gives the filters zero prehistory, even for very short paths.
  noise <- c(rep(0, length(ar) + length(ma)), start.innov, innov)
  if (length(ma)) {
    noise <- stats::filter(noise, c(1, ma), sides = 1)
    noise[seq_along(ma)] <- 0
  }
  if (length(ar)) noise <- stats::filter(noise, ar, method = "recursive")
  series <- as.numeric(utils::tail(noise, nsim))
  d <- object$arma[6L]
  D <- object$arma[7L]
  m <- object$arma[5L]
  if (d > 0L) series <- utils::tail(stats::diffinv(series, differences = d), nsim)
  if (D > 0L) series <- utils::tail(stats::diffinv(series, lag = m, differences = D), nsim)
  if ("intercept" %in% names(object$coef)) series <- series + object$coef["intercept"]
  if ("drift" %in% names(object$coef)) series <- series + object$coef["drift"] * seq_len(nsim)
  stats::ts(series, frequency = object$frequency)
}

arima_check_stationarity <- function(object) {
  decay <- Inf
  for (prefix in c("ar", "sar")) {
    coefficients <- object$coef[grep(paste0("^", prefix, "[0-9]+$"), names(object$coef))]
    if (any(!is.finite(coefficients))) stop("AR coefficients must be finite")
    if (!any(coefficients != 0)) next
    roots <- polyroot(c(1, -coefficients))
    if (any(Mod(roots) <= 1)) stop("AR part of model must be stationary")
    period <- if (prefix == "sar") object$frequency else 1
    decay <- min(decay, log(min(Mod(roots))) / period)
  }
  decay
}

arima_burnin <- function(object) {
  decay <- arima_check_stationarity(object)
  ar <- head(object$model$phi, max(c(0L, which(object$model$phi != 0))))
  ma <- head(object$model$theta, max(c(0L, which(object$model$theta != 0))))
  orders <- length(ar) + length(ma)
  if (!length(ar)) return(orders)
  size <- 2^ceiling(log2(max(64, 2 * (orders + 1), ceiling(12 / decay))))
  # Normalize the impulse response before squaring to avoid energy overflow.
  while (size <= 2^22) {
    response <- c(1, ma, rep(0, size - length(ma) - 1L))
    response <- as.numeric(stats::filter(response, ar, method = "recursive"))
    if (any(!is.finite(response))) break
    energy <- (response / max(abs(response)))^2
    if (sum(energy[seq.int(size / 2 + 1L, size)]) <= 1e-8 * sum(energy)) {
      return(size / 2)
    }
    size <- size * 2
  }
  stop("automatic burn-in did not converge; supply n.start explicitly")
}
