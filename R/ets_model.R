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
#'
#' @param frequency The length of the seasonal period (e.g., 12 for monthly data).
#' @param error Character string specifying the error component: \code{"A"} for
#'   additive or \code{"M"} for multiplicative.
#' @param trend Character string specifying the trend component: \code{"N"} for
#'   none or \code{"A"} for additive.
#' @param seasonal Character string specifying the seasonal component:
#'   \code{"N"} for none, \code{"A"} for additive, or \code{"M"} for
#'   multiplicative.
#' @param damped Logical value indicating whether an additive trend is damped.
#' @param alpha Smoothing parameter for the level.
#' @param beta Smoothing parameter for the trend.
#' @param gamma Smoothing parameter for seasonality.
#' @param phi Damping parameter.
#' @param level Initial level \eqn{\ell_0}.
#' @param slope Initial slope \eqn{b_0}.
#' @param season Numeric vector of initial seasonal states
#'   \eqn{s_{1-m},\dots,s_0}.
#' @param sigma The standard deviation of the noise.
#' @return An \code{ets} object compatible with the \pkg{forecast} package's
#'   simulation methods.
#' @author Rob J Hyndman
#' @seealso \code{\link[forecast]{simulate.ets}}
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
  # error
  if (is.null(error)) {
    error <- sample(c("A", "M"), 1)
  }
  # sigma
  if (is.null(sigma)) {
    if (error == "A") { # Choose variance in [1,5]
      sigma2 <- runif(1, 1, 5)
    } else {
      sigma2 <- runif(1, 1e-4, 0.05)
    }
  } else {
    sigma2 <- sigma^2
  }
  # trend
  if (is.null(trend)) {
    trend <- sample(c("N", "A"), 1)
  }
  if (is.null(damped)) {
    if (trend == "A") {
      damped <- sample(c(TRUE, FALSE), 1)
    } else {
      damped <- FALSE
    }
  }
  # seasonal
  if (is.null(seasonal)) {
    if (frequency > 1) {
      seasonal <- sample(c("N", "A", "M"), 1)
    } else {
      seasonal <- "N"
    }
  }
  if (seasonal != "N") {
    if (frequency <= 1) {
      stop("Seasonal models must have frequency greater than 1")
    } else {
      m <- frequency
    }
  } else {
    m <- 1
  }

  # Find admissible parameters
  admissible <- FALSE
  while (!admissible) {
    # alpha
    if (is.null(alpha)) {
      alpha <- runif(1, 0, 1)
    }
    # beta
    if (trend != "N" & is.null(beta)) {
      beta <- runif(1, 0, alpha)
    }
    # gamma
    if (seasonal != "N" & is.null(gamma)) {
      gamma <- runif(1, 0, 1 - alpha)
    }
    # phi
    if (trend == "A" & damped & is.null(phi)) {
      phi <- runif(1, 0.8, 0.98)
    }
    admissible <- ets_admissible(alpha, beta, gamma, phi, m)
  }
  # Find initial states
  if (is.null(level)) {
    if (error == "A") {
      level <- runif(1, -1, 1)
    } else {
      level <- runif(1, 2, 10)
    }
  }
  if (trend != "N" & is.null(slope)) {
    slope <- runif(1, -1, 1)
  }
  if (is.null(season)) {
    if (seasonal == "A") {
      season <- runif(m - 1, -1, 1)
      season <- c(season, -sum(season))
    } else if (seasonal == "M") {
      season <- runif(m - 1, 0.5, 1.5)
      season <- c(season, m - sum(season))
    }
    if (seasonal != "N") {
      season <- matrix(season, nrow = 1)
      colnames(season) <- paste0("s", seq(m))
    }
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
  if (phi < 0 || phi > 1 + 1e-8) {
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
    roots <- polyroot(P)

    if (max(abs(roots)) > 1 + 1e-10) {
      return(FALSE)
    }
  }
  # Passed all tests
  return(TRUE)
}
