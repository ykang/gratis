#' Generate random MAR time series
#'
#' This function is deprecated. Use \code{\link{mar_model}()} with
#' \code{\link{generate.mar}()} or \code{\link{simulate.mar}()} instead.
#'
#' @param n.ts Number of time series to generate.
#' @param freq Seasonal period of the generated series.
#' @param nComp Number of mixture components. If \code{NULL}, this is sampled
#'   randomly for each generated series.
#' @param n Length of each generated series.
#' @param output_format Output format: \code{"list"} for the historical list
#'   format or \code{"tsibble"} for a list of tsibbles.
#' @return When \code{output_format = "list"}, a list containing generated time
#'   series and the parameters used for each series. When
#'   \code{output_format = "tsibble"}, a list of tsibbles.
#' @author Yanfei Kang and Feng Li
#' @examples
#' x <- generate_ts(n.ts = 2, freq = 12, nComp = 2, n = 120)
#' x$N1$pars
#' forecast::autoplot(x$N1$x)
#' @references Wong, CS & WK Li (2000).
#' @export
generate_ts <- function(n.ts = 1, freq = 1, nComp = NULL, n = 120, output_format = "list") {
  warning("This function is deprecated. It is recommended you use mar_model() and simulate.mar() or generate.mar() instead.")
  count <- 1
  generated.mixture.data <- list()
  sigmas <- sample(c(1:5), 5, replace = TRUE)
  phi0 <- sample(c(0, 3), 1)
  while (count <= n.ts) {
    if (is.null(nComp)) {
      # set the number of components
      nComp <- sample(1:5, 1)
    } else {
      nComp <- nComp
    }
    pars <- list()
    # seasonal data
    if (freq != 1) {
      for (i in 1:nComp) {
        pars[[sprintf("pars%d", i)]] <- rnorm(4, 0, 0.5)
      }
      mixture.weights <- rep(NA, nComp)
      for (i in 1:nComp) {
        mixture.weights[i] <- runif(1)
      }
      mixture.weights <- mixture.weights / sum(mixture.weights)
      means.ar.par.list <- lapply(pars, function(x) {
        d <- sample(c(0, 1), 1, prob = c(0.1, 0.9))
        D <- sample(c(0, 1), 1, prob = c(0.6, 0.4))
        c(phi0, pi_coefficients(ar = x[1:2], sar = x[3:4], d = d, D = D, m = freq))
      })
      sigmas.list <- list()
      for (i in 1:nComp) {
        sigmas.list[[i]] <- rep(sigmas[i], n + freq * 10)
      }
      pars$weights <- mixture.weights
      x <- rmixnorm_ts(
        n = n + freq * 10,
        means.ar.par.list = means.ar.par.list,
        sigmas.list = sigmas.list,
        weights = mixture.weights
      )
      # allow spikes
      if (runif(1) <= 0.01) {
        x[sample(1:length(x), 1)] <- max(x) * sample(2:5, 1)
      }
      x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
      if (!(any(is.na(x)) || (max(abs(x), na.rm = TRUE) > 1e5))) {
        number <- paste0("N", count)
        generated.mixture.data[[number]] <- list()
        generated.mixture.data[[number]]$x <- x
        generated.mixture.data[[number]]$pars <- pars
        count <- count + 1
      }
    } else {
      # nonseasonal data
      for (i in 1:nComp) {
        pars[[sprintf("pars%d", i)]] <- rnorm(2, 0, 0.5)
      }
      mixture.weights <- rep(NA, nComp)
      for (i in 1:nComp) {
        mixture.weights[i] <- runif(1)
      }
      mixture.weights <- mixture.weights / sum(mixture.weights)
      means.ar.par.list <- lapply(pars, function(x) {
        d <- sample(c(0, 1, 2), 1, prob = c(0.1, 0.6, 0.3))
        c(phi0, pi_coefficients(ar = x[1:2], d = d, m = freq))
      })
      sigmas.list <- list()
      for (i in 1:nComp) {
        sigmas.list[[i]] <- rep(sigmas[i], n + freq * 10)
      }
      pars$weights <- mixture.weights
      x <- rmixnorm_ts(
        n = n + freq * 10,
        means.ar.par.list = means.ar.par.list,
        sigmas.list = sigmas.list,
        weights = mixture.weights
      )
      # allow spikes
      if (runif(1) <= 0.01) {
        x[sample(1:length(x), 1)] <- max(x) * sample(2:5, 1)
      }
      x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
      if (max(abs(x), na.rm = TRUE) < 1e5) {
        number <- paste0("N", count)
        generated.mixture.data[[number]] <- list()
        generated.mixture.data[[number]]$x <- x
        generated.mixture.data[[number]]$pars <- pars
        count <- count + 1
      }
    }
  }
  # New content
  output <- if (output_format == "list") {
    generated.mixture.data
  } else if (output_format == "tsibble") {
    x <- generated.mixture.data
    map(x, ~ as_tsibble(.x$x))
  }
  return(output)
}



#' Generate a random multiple-seasonal MAR time series
#'
#' This function is deprecated. Use \code{\link{mar_model}()} with
#' \code{\link{generate.mar}()} or \code{\link{simulate.mar}()} instead.
#'
#' @param seasonal.periods Numeric vector of seasonal periods.
#' @param n Length of the generated series.
#' @param nComp Number of mixture components. If \code{NULL}, this is sampled
#'   randomly for each generated component series.
#' @param output_format Output format: \code{"list"} for the historical
#'   \code{msts} result or \code{"tsibble"} for a tsibble.
#' @return A multiple-seasonal time series when \code{output_format = "list"},
#'   or a tsibble when \code{output_format = "tsibble"}.
#' @export
#' @examples
#' x <- generate_msts(seasonal.periods = c(7, 365), n = 800, nComp = 2, output_format = "list")
#' forecast::autoplot(x)
generate_msts <- function(seasonal.periods = c(7, 365), n = 800, nComp = NULL, output_format = "list") {
  warning("This function is deprecated. It is recommended you use mar_model() and simulate.mar() or generate.mar() instead.")
  x.list <- map(seasonal.periods, function(p) {
    suppressWarnings(generate_ts(n.ts = 1, freq = p, n = n, nComp = nComp)$N1$x)
  })
  names(x.list) <- paste0("Season", seasonal.periods)
  x.list[1:(length(x.list) - 1)] <- lapply(x.list[1:(length(x.list) - 1)], function(x) {
    x - trendcycle(stl(x, "per"))
  })
  weights <- msts_weights(length(seasonal.periods))
  res <- as_tibble(scale(x.list %>% bind_cols())[, ]) %>%
    mapply("*", ., weights) %>%
    as_tibble() %>%
    mutate(x = rowSums(.)) %>%
    select(x) %>%
    msts(seasonal.periods = seasonal.periods)
  # New content
  output <- if (output_format == "list") {
    res
  } else if (output_format == "tsibble") {
    tsibble::as_tsibble(res)
  }
  return(output)
}

# ===========================================================================
# Simulated weights for the simulation of msts
# ===========================================================================
msts_weights <- function(n.periods) {
  gamma <- runif(n.periods, 0)
  weights <- gamma / sum(gamma)
  return(weights)
}
utils::globalVariables(c(".", "x", "Season2"))
