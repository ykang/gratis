#' Specify an ETS model
#'
#' Construct an ETS state space model that can be used with
#' \code{\link[forecast]{simulate.ets}()} or \code{\link{generate.ets}()}.
#' Any omitted argument is randomly selected from the supported parameter space.
#'
#' The error component is sampled from \code{"A"} and \code{"M"}, the trend
#' component from \code{"N"} and \code{"A"}, and the seasonal component from
#' \code{"N"}, \code{"A"}, and \code{"M"} when \code{frequency > 1}. Smoothing
#' parameters are sampled until they satisfy the ETS admissibility conditions.
#' The innovation variance is sampled from \eqn{(1,5)} for additive errors and
#' from \eqn{(0.0001,0.05)} for multiplicative errors. Initial states are sampled
#' from ranges appropriate to the selected components.
#' Omitted smoothing parameters are resampled after rejection, with a maximum
#' of 10000 attempts. Fixed parameters must satisfy the usual smoothing bounds:
#' \eqn{0 < \alpha < 1}, \eqn{0 \le \beta \le \alpha},
#' \eqn{0 \le \gamma \le 1-\alpha}, and \eqn{0 < \phi \le 1}.
#' Multiplicative seasonal factors are positive and sampled factors are
#' normalized to have mean one.
#'
#' Simulation starts from the supplied initial states without an automatic
#' burn-in. Multiplicative innovations are clipped below at -1 by the
#' \pkg{forecast} simulator, matching Python's convention. Admissibility and
#' positive initial states do not guarantee positive or finite values for every
#' simulated path.
#'
#' @param frequency Positive integer seasonal period, for example 12 for monthly
#'   data. Seasonal models require a value greater than 1.
#' @param error Character string specifying the error component: \code{"A"} for
#'   additive or \code{"M"} for multiplicative.
#' @param trend Character string specifying the trend component: \code{"N"} for
#'   none or \code{"A"} for additive.
#' @param seasonal Character string specifying the seasonal component:
#'   \code{"N"} for none, \code{"A"} for additive, or \code{"M"} for
#'   multiplicative.
#' @param damped Logical value indicating whether an additive trend is damped.
#' @param alpha Level smoothing parameter, strictly between 0 and 1.
#' @param beta Trend smoothing parameter in \eqn{[0,\alpha]}. Ignored without a trend.
#' @param gamma Seasonal smoothing parameter in \eqn{[0,1-\alpha]}. Ignored without seasonality.
#' @param phi Damping parameter in \eqn{(0,1]}. Ignored without a damped trend.
#' @param level Finite initial level \eqn{\ell_0}. Must be positive for
#'   multiplicative errors or seasonality.
#' @param slope Finite initial additive trend state \eqn{b_0}. Ignored without a trend.
#' @param season Finite numeric vector of initial seasonal states
#'   \eqn{s_{1-m},\dots,s_0}, ordered oldest to newest, of length \code{frequency}.
#'   Multiplicative factors must be positive. Supplied states are not normalized.
#' @param sigma The standard deviation of the noise. Its squared variance must
#'   lie in the positive normal floating-point range, in both R and Python.
#'   It is measured on the observation scale for additive errors and is
#'   dimensionless for multiplicative errors.
#' @return An \code{ets} object compatible with the \pkg{forecast} package's
#'   simulation methods.
#' @author Rob J Hyndman
#' @seealso \code{\link[forecast]{simulate.ets}}, \code{\link{generate.ets}}
#' @examples
#' # An ETS(A,A,N) model with random parameters
#' model1 <- ets_model(error = "A", trend = "A", seasonal = "N")
#' # An ETS(A,A,N) model with specific parameters
#' model2 <- ets_model(
#'   error = "A", trend = "A", seasonal = "N",
#'   alpha = 0.3, beta = 0.2, level = 0, slope = 1, sigma = 2
#' )
#' # A quarterly seasonal ETS model with random admissible parameters
#' model3 <- ets_model(error = "A", trend = "N", seasonal = "A", frequency = 4)
#' # Simulate from each model and plot the results
#' plot(simulate(model1, 100))
#' plot(simulate(model2, 100))
#' plot(simulate(model3, 100))
#' @export
ets_model <- function(frequency = 1, error = NULL, trend = NULL, seasonal = NULL,
                      alpha = NULL, beta = NULL, gamma = NULL, phi = NULL,
                      level = NULL, slope = NULL, season = NULL, damped = NULL,
                      sigma = NULL) {
  if (length(frequency) != 1L || !is.finite(frequency) ||
      frequency < 1 || frequency != as.integer(frequency)) {
    stop("frequency must be a positive integer")
  }
  # error
  if (is.null(error)) {
    error <- sample(c("A", "M"), 1)
  }
  error <- match.arg(toupper(error), c("A", "M"))
  # sigma
  if (is.null(sigma)) {
    if (error == "A") { # Choose variance in [1,5]
      sigma2 <- runif(1, 1, 5)
    } else {
      sigma2 <- runif(1, 1e-4, 0.05)
    }
  } else {
    validate_sigma(sigma)
    sigma2 <- sigma^2
  }
  # trend
  if (is.null(trend)) {
    trend <- sample(c("N", "A"), 1)
  }
  trend <- match.arg(toupper(trend), c("N", "A"))
  if (is.null(damped)) {
    if (trend == "A") {
      damped <- sample(c(TRUE, FALSE), 1)
    } else {
      damped <- FALSE
    }
  }
  if (trend == "N") damped <- FALSE
  # seasonal
  if (is.null(seasonal)) {
    if (frequency > 1) {
      seasonal <- sample(c("N", "A", "M"), 1)
    } else {
      seasonal <- "N"
    }
  }
  seasonal <- match.arg(toupper(seasonal), c("N", "A", "M"))
  if (seasonal != "N") {
    if (frequency <= 1) {
      stop("Seasonal models must have frequency greater than 1")
    } else {
      m <- frequency
    }
  } else {
    m <- 1
  }

  if (trend == "N") beta <- NULL
  if (seasonal == "N") gamma <- NULL
  if (!damped) phi <- NULL
  supplied <- list(alpha = alpha, beta = beta, gamma = gamma, phi = phi)
  for (name in names(supplied)) {
    value <- supplied[[name]]
    if (is.null(value)) next
    if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1) {
      stop(paste(name, "must be finite and between 0 and 1"))
    }
    if (name %in% c("alpha", "phi") && value == 0) {
      stop(paste(name, "must be positive"))
    }
    if (name == "alpha" && value == 1) stop("alpha must be less than 1")
  }
  lower <- if (is.null(beta)) 0 else beta
  upper <- if (is.null(gamma)) 1 else 1 - gamma
  if (lower > upper || upper <= 0 || lower >= 1) {
    stop("beta and gamma leave no admissible alpha")
  }
  if (!is.null(alpha) && (alpha < lower || alpha > upper)) {
    stop("beta must be <= alpha and gamma must be <= 1 - alpha")
  }
  all_fixed <- !is.null(alpha) && (trend == "N" || !is.null(beta)) &&
    (seasonal == "N" || !is.null(gamma)) && (!damped || !is.null(phi))
  admissible <- FALSE
  for (attempt in seq_len(10000L)) {
    a <- if (is.null(alpha)) runif(1, lower, upper) else alpha
    b <- if (trend == "N") NULL else if (is.null(beta)) runif(1, 0, a) else beta
    g <- if (seasonal == "N") NULL else if (is.null(gamma)) runif(1, 0, 1 - a) else gamma
    f <- if (!damped) NULL else if (is.null(phi)) runif(1, 0.8, 0.98) else phi
    admissible <- ets_admissible(a, b, g, f, m)
    if (admissible) break
    if (all_fixed) stop("inadmissible ETS smoothing parameters")
  }
  if (!admissible) stop("could not find admissible ETS smoothing parameters")
  alpha <- a
  beta <- b
  gamma <- g
  phi <- f
  # Find initial states
  if (is.null(level)) {
    if (error == "A" && seasonal != "M") {
      level <- runif(1, -1, 1)
    } else {
      level <- runif(1, 2, 10)
    }
  }
  if (length(level) != 1L || !is.finite(level)) stop("level must be finite")
  if ((error == "M" || seasonal == "M") && level <= 0) {
    stop("multiplicative ETS models require a positive initial level")
  }
  if (trend != "N" & is.null(slope)) {
    slope <- runif(1, -1, 1)
  }
  if (trend == "N") slope <- NULL
  if (!is.null(slope) && (length(slope) != 1L || !is.finite(slope))) {
    stop("slope must be finite")
  }
  if (is.null(season)) {
    if (seasonal == "A") {
      season <- runif(m - 1, -1, 1)
      season <- c(season, -sum(season))
    } else if (seasonal == "M") {
      season <- runif(m, 0.5, 1.5)
      season <- season / mean(season)
    }
  }
  if (seasonal != "N") {
    if (length(season) != m || any(!is.finite(season))) {
      stop("season must be finite and have length equal to frequency")
    }
    if (seasonal == "M" && any(season <= 0)) {
      stop("multiplicative seasonal factors must be positive")
    }
    # forecast stores seasonal states from newest to oldest.
    season <- matrix(rev(season), nrow = 1)
    colnames(season) <- paste0("s", seq(m))
  } else {
    season <- NULL
  }
  states <- cbind(l = level, b = slope, season)
  # Components vector
  components <- c(error, trend, seasonal, as.character(damped))
  # Method string
  method <- paste0("ETS(", error, ",", trend, ifelse(damped, "d", ""), ",", seasonal, ")")
  # Return
  structure(list(
    method = method,
    states = states, initstate = states, components = components, m = frequency,
    par = c(alpha = alpha, beta = beta, gamma = gamma, phi = phi), sigma2 = sigma2
  ), class = "ets")
}

# Check admissibility (borrowed from forecast package)

ets_admissible <- function(alpha, beta, gamma, phi, m) {
  if (is.null(phi)) {
    phi <- 1
  }
  if (any(!is.finite(c(alpha, beta, gamma, phi))) || phi <= 0 || phi > 1 + 1e-8) {
    return(FALSE)
  }
  if (is.null(gamma)) {
    if (alpha < 1 - 1 / phi || alpha > 1 + 1 / phi) {
      return(FALSE)
    }
    if (!is.null(beta)) {
      if (beta < alpha * (phi - 1) || beta > (1 + phi) * (2 - alpha)) {
        return(FALSE)
      }
    }
  } else if (m > 1) { # Seasonal model
    if (is.null(beta)) {
      beta <- 0
    }
    if (gamma < max(1 - 1 / phi - alpha, 0) || gamma > 1 + 1 / phi - alpha) {
      return(FALSE)
    }
    if (alpha < 1 - 1 / phi - gamma * (1 - m + phi + phi * m) / (2 * phi * m)) {
      return(FALSE)
    }
    if (beta < -(1 - phi) * (gamma / m + alpha)) {
      return(FALSE)
    }

    # End of easy tests. Now use characteristic equation
    P <- c(
      phi * (1 - alpha - gamma), alpha + beta - alpha * phi + gamma - 1,
      rep(alpha + beta - alpha * phi, m - 2), (alpha + beta - phi), 1
    )
    # Use the companion eigenvalues, as in NumPy, for long seasonal periods.
    companion <- rbind(-rev(P[-length(P)]), cbind(diag(length(P) - 2L), 0))
    roots <- eigen(companion, only.values = TRUE)$values

    if (max(abs(roots)) > 1 + 1e-10) {
      return(FALSE)
    }
  }
  # Passed all tests
  return(TRUE)
}
