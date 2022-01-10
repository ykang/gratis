# Fitness function for time series generation.
#
# @param pars Parameters
# @param features Time series features.
# @param seasonal Seasonal effects.
# @param n Length of time series
# @param freq Frequence of time series
# @param target Target time series features
# @param nComp No. of components used in mixture models.
# @param selected.features Selected features.
#
# @return NA
#' @importFrom forecast ndiffs

fitness_ts <- function(pars, features, seasonal, n = 120, freq = 12, target, nComp, selected.features) {
  pars <- pars2list(pars, seasonal, nComp)
  if (seasonal < 2) {
    x <- pars2x(pars, seasonal, freq, nComp, n)
    if (max(abs(x)) > 1e5) {
      return(list(value = -100, x = x))
    }
    return(list(
      # value = -sqrt(sum((tsfeatures(x, features = features) %>% select(selected.features) - target)^2)) / sqrt(sum(target^2)),
      value = -sqrt(sum((tsfeatures(x, features = features) %>% select(selected.features) - target)^2)),
      x = x
    ))
  } else {
    x.list <- as.list(rep(0, length(freq)))
    for (i in seq(freq)) {
      x.list[[i]] <- pars2x(pars[[i]], 1, freq[i], nComp, n)
    }
    x.list[1:(length(x.list) - 1)] <- lapply(x.list[1:(length(x.list) - 1)], function(x) {
      x - trendcycle(stl(x, "per"))
    })
    names(x.list) <- paste0("Season", seq(freq))
    res <- as_tibble(scale(x.list %>% bind_cols())[, ]) %>% mutate(Season2 = Season2 / pars$ms.scale)
    x <- res %>%
      mutate(x = rowSums(.[seq(freq)])) %>%
      dplyr::select(x) %>%
      unlist() %>%
      as.numeric() %>%
      msts(seasonal.periods = freq)
    if (is.na(max(abs(x))) | max(abs(x)) > 1e5) {
      return(list(value = -100, x = x))
    }
    return(list(
      # value = -sqrt(sum((tsfeatures::tsfeatures(x, features = features) %>%
      #   select(selected.features) - target)^2)) / sqrt(sum(target^2)),
      value = -sqrt(sum((tsfeatures::tsfeatures(x, features = features) %>%
        select(selected.features) - target)^2)),
      x = x
    ))
  }
}

pars2x <- function(pars, seasonal, freq, nComp, n) {
  if (seasonal == 0) {
    means.ar.par.list <- lapply(pars[1:nComp], function(x) {
      c(0, pi_coefficients(ar = x[1:2], m = freq))
    })
    sigmas.list <- list()
    for (i in 1:nComp) {
      sigmas.list[[i]] <- rep(i, n + freq * 10)
    }
    weights <- pars$weights
    x <- rmixnorm_ts(
      n = n + freq * 10,
      means.ar.par.list = means.ar.par.list,
      sigmas.list = sigmas.list,
      weights = weights
    )
    x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
    if (pars$p.change <= 0.1) {
      t.change <- sample(1:n, 1)
      x[t.change:n] <- x[t.change:n] + 50
    }
    if (pars$p.diff <= 0.5) {
      x <- ts(diffinv(x), frequency = freq)
      x <- window(x, start = c(1, 2))
    }
  } else {
    means.ar.par.list <- lapply(pars[1:nComp], function(x) {
      # d = sample(c(0,1), 1)
      # D =  sample(c(0,1), 1)
      c(0, pi_coefficients(ar = x[1:2], sar = x[3:4], d = 0, D = 0, m = freq))
    })
    sigmas.list <- list()
    for (i in 1:nComp) {
      sigmas.list[[i]] <- rep(i, n + freq * 10)
    }
    weights <- pars$weights
    x <- rmixnorm_ts(
      n = n + freq * 10,
      means.ar.par.list = means.ar.par.list,
      sigmas.list = sigmas.list,
      weights = weights
    )
    x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
    if (pars$p.change <= 0.1) {
      t.change <- sample(1:n, 1)
      x[t.change:n] <- x[t.change:n] + 50
    }
    if (pars$p.diff <= 0.5) {
      x <- ts(diffinv(x), frequency = freq)
      x <- window(x, start = c(1, 2))
    }
    if (pars$p.Diff <= 0.5) {
      x <- ts(diffinv(x, lag = freq), frequency = freq)
      x <- window(x, start = c(2, 1))
    }
  }
  return(x)
}

pars2x1 <- function(pars, seasonal, freq, nComp, n, x0) {
  if (seasonal == 0) {
    means.ar.par.list <- lapply(pars[1:nComp], function(x) {
      c(0, pi_coefficients(ar = x[1:2], m = freq, tol = 1e-07))
    })
    sigmas.list <- list()
    for (i in 1:nComp) {
      sigmas.list[[i]] <- rep(i, n + freq * 10)
    }
    weights <- pars$weights
    x <- rmixnorm_ts(
      n = n + freq * 10,
      means.ar.par.list = means.ar.par.list,
      sigmas.list = sigmas.list,
      weights = weights
    )
    x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
    if (pars$p.change <= 0) {
      t.change <- sample(1:n, 1)
      x[t.change:n] <- x[t.change:n] + 50
    }
    if (pars$p.diff <= ifelse(ndiffs(x0) > 0, 0.5, 0)) {
      x <- ts(diffinv(x, differences = ndiffs(x0)), frequency = freq)
      x <- window(x, start = c(1, 1 + ndiffs(x0)))
    }
  } else {
    means.ar.par.list <- lapply(pars[1:nComp], function(x) {
      # d = sample(c(0,1), 1)
      # D =  sample(c(0,1), 1)
      c(0, pi_coefficients(ar = x[1:2], sar = x[3:4], d = 0, D = 0, m = freq))
    })
    sigmas.list <- list()
    for (i in 1:nComp) {
      sigmas.list[[i]] <- rep(i, n + freq * 10)
    }
    weights <- pars$weights
    x <- rmixnorm_ts(
      n = n + freq * 10,
      means.ar.par.list = means.ar.par.list,
      sigmas.list = sigmas.list,
      weights = weights
    )
    x <- ts(x[(1 + freq * 10):(n + freq * 10)], frequency = freq)
    if (pars$p.change <= 0) {
      t.change <- sample(1:n, 1)
      x[t.change:n] <- x[t.change:n] + 50
    }
    if (pars$p.diff <= ifelse(forecast::ndiffs(x0) > 0, 0.25, 0)) {
      x <- ts(diffinv(x), frequency = freq)
      x <- window(x, start = c(1, 2))
    }
    if (pars$p.Diff <= ifelse(forecast::nsdiffs(x0) > 0, 0.25, 0)) {
      x <- ts(diffinv(x, lag = freq), frequency = freq)
      x <- window(x, start = c(2, 1))
    }
  }
  return(x)
}


pars2list <- function(pars, seasonal, nComp) {
  parslist <- list()
  if (seasonal == 0) {
    for (i in 1:nComp) {
      parslist[[sprintf("pars%d", i)]] <- pars[(2 * i - 1):(2 * i)]
    }
    parslist$p.change <- pars[2 * nComp + 1]
    parslist$p.diff <- pars[2 * nComp + 2]
    parslist$weights <- c(pars[(2 * nComp + 3):(3 * nComp + 1)], 1 - sum(pars[(2 * nComp + 3):(3 * nComp + 1)]))
    parslist$weights <- exp(parslist$weights) / sum(exp(parslist$weights))
  } else if (seasonal == 1) {
    for (i in 1:nComp) {
      parslist[[sprintf("pars%d", i)]] <- pars[(4 * i - 3):(4 * i)]
    }
    parslist$p.change <- pars[4 * nComp + 1]
    parslist$p.diff <- pars[4 * nComp + 2]
    parslist$p.Diff <- pars[4 * nComp + 3]
    parslist$weights <- c(pars[(4 * nComp + 4):(5 * nComp + 2)], 1 - sum(pars[(4 * nComp + 4):(5 * nComp + 2)]))
    parslist$weights <- exp(parslist$weights) / sum(exp(parslist$weights))
  } else {
    pars1 <- pars[1:((length(pars) - 1) / 2)]
    pars2 <- pars[((length(pars) - 1) / 2 + 1):(length(pars) - 1)]
    parslist[[1]] <- pars2list(pars1, 1, nComp)
    parslist[[2]] <- pars2list(pars2, 1, nComp)
    parslist$ms.scale <- pars[length(pars)]
  }
  return(parslist)
}

#  scalets(), scale time series function is adapted from tsfeatures package:
#  tsfeatures: Time Series Feature Extraction. R package version 1.0.2.
#  https://CRAN.R-project.org/package=tsfeatures
scalets <- function(x) {
  n <- length(x)
  if (forecast::is.constant(x)) {
    return(x)
  }
  scaledx <- as.numeric(scale(x, center = TRUE, scale = TRUE))
  if ("msts" %in% class(x)) {
    msts <- attributes(x)$msts
    y <- forecast::msts(scaledx, seasonal.periods = msts)
  } else {
    y <- as.ts(scaledx)
  }
  tsp(y) <- tsp(x)
  return(y)
}

fitness_ts1 <- function(pars, x0, seasonal, n = 60, freq = 12, nComp, h = 18) {
  pars <- pars2list(pars, seasonal, nComp)
  if (seasonal < 2) {
    x <- pars2x1(pars, seasonal, freq, nComp, n, x0)
    x1 <- x[1:(length(x) - h)]
    x2 <- x[(length(x) - h + 1):length(x)]
    x <- scalets(x1)
    xx <- (x2 - mean(x1)) / sd(x1)
    if (max(abs(x)) > 1e5) {
      return(list(value = -100, x = c(x, xx)))
    }
    return(list(
      # pearson correlation distance
      # value = cor(as.vector(x[1:length(x0)]), as.vector(x0)),
      # cort distance
      # value = - diss.cort(as.vector(x), as.vector(x0), k = 2),
      # value = tsgeneration:::corrtemporder1(as.vector(x), as.vector(x0)),
      # value = - mean(abs(as.vector(x[1:length(x0)]) - as.vector(x0))),
      value = -sqrt(sum((as.vector(x) - as.vector(x0))^2)),
      x = c(x, xx)
    ))
  } else {
    x.list <- as.list(rep(0, length(freq)))
    for (i in seq(freq)) {
      x.list[[i]] <- pars2x1(pars[[i]], 1, freq[i], nComp, n)
    }
    x.list[1:(length(x.list) - 1)] <- lapply(x.list[1:(length(x.list) - 1)], function(x) {
      x - trendcycle(stl(x, "per"))
    })
    names(x.list) <- paste0("Season", seq(freq))
    res <- as_tibble(scale(x.list %>% bind_cols())[, ]) %>% mutate(Season2 = Season2 / pars$ms.scale)
    x <- res %>%
      mutate(x = rowSums(.[seq(freq)])) %>%
      dplyr::select(x) %>%
      unlist() %>%
      as.numeric() %>%
      msts(seasonal.periods = freq)
    if (is.na(max(abs(x))) | max(abs(x)) > 1e5) {
      return(list(value = -100, x = x))
    }
    return(list(
      # value = -sqrt(sum((tsfeatures::tsfeatures(x, features = features) %>%
      #   select(selected.features) - target)^2)) / sqrt(sum(target^2)),
      value = -sqrt(sum((scalets(as.vector(x)) -
        scalets(as.vector(x0)))^2)),
      x = x
    ))
  }
}
