#' Specify parameters for an ETS model
#'
#' This function allows the parameters of a ETS state space model to be specified.
#' The output can be used in \code{\link[forecast]{simulate.ets}()}
#' and \code{\link{generate.ets}}.
#' If any argument is \code{NULL}, the corresponding parameters are randomly selected.
#' The error component is chosen from \{A,M\}, the trend component is chosen from
#' \{N,A,Ad\}, and the seasonal component is chosen from \{N,A,M\}. In all cases,
#' the component is selected uniformly on the options. The parameters are selected
#' uniformly on the forecastable parameter space. The noise variance sigma
#' is uniformly sampled on (1,5) for additive errors, and on (0.0001,0.05) for
#' multiplicative errors. The initial states are chosen uniformly on (-1,1)
#' in all cases except for multiplicative seasonal states which are uniform on
#' (0.5, 1.5), and models with multiplicative errors for which the level is uniform
#' on (2, 10). The parameterization is as specified in Hyndman & Athanasopoulos (2021).
#'
#' @param frequency The length of the seasonal period (e.g., 12 for monthly data).
#' @param error A character string specifying the error part of the ETS model: either "A" or "M".
#' @param trend A character string specifying the trend part of the ETS model: either "N", "A" or "Ad".
#' @param seasonal A character string specifying the seasonal part of the ETS model: either "N", "A" or "M".
#' @param damped A logical value indicating if the trend is damped or not.
#' @param alpha A numeric value for the smoothing parameter controlling the level.
#' @param beta A numeric value for the smoothing parameter controlling the trend.
#' @param gamma A numeric value for the smoothing parameter controlling the seasonality.
#' @param phi A numeric value specifying the damping parameter.
#' @param level A numeric value specifying the initial level \eqn{\ell_0}.
#' @param slope A numeric value specifying the initial slope \eqn{b_0}
#' @param season A numeric vector specifying the initial states \eqn{s_{1-m},...,s_0}.
#' @param sigma The standard deviation of the noise.
#' @return An `ets` object as described in the \code{\link[forecast]{ets}} function from the forecast package.
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
#' # A multiplicative quarterly seasonal ETS model with random parameters
#' model3 <- ets_model(seasonal = "M", frequency = 4)
#' # Simulate from each model and plot the results
#' library(forecast)
#' simulate(model1, 100) %>% plot()
#' simulate(model2, 100) %>% plot()
#' simulate(model3, 100) %>% plot()
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
