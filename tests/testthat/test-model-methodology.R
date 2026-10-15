test_that("AR expansion uses the positive MA convention and keeps long lags", {
  expect_equal(head(pi_coefficients(ma = 0.5), 6), 0.5 * (-0.5)^(0:5))
  expect_equal(pi_coefficients(ar = 0.8, ma = -0.8), numeric(0))
  expected <- numeric(731)
  expected[c(1, 365, 366, 730, 731)] <- c(0.3, 0.2, -0.06, 0.3, -0.09)
  expect_equal(pi_coefficients(ar = 0.3, sar = c(0.2, 0.3), m = 365), expected)
  model <- mar_model(phi = matrix(0.3), Phi = matrix(c(0.2, 0.3)), d = 0,
                     D = 0, constants = 0, sigmas = 1, weights = 1,
                     seasonal_periods = 365)
  expect_equal(as.numeric(model$ar), c(0, expected))
  expect_warning(pi_coefficients(ma = 0.99, n_terms = 10), "increase n_terms")
})

test_that("ARIMA constant represents a mean or linear drift", {
  model <- arima_model(p = 1, d = 0, q = 0, phi = 0.5, constant = 2, sigma = 1e-8)
  expect_equal(as.numeric(model$coef), c(0.5, 2))
  baseline <- arima_model(p = 1, d = 0, q = 0, phi = 0.5, constant = 0, sigma = 1e-8)
  set.seed(1)
  x <- as.numeric(simulate(model, nsim = 10, innov = rep(0, 10)))
  set.seed(1)
  y <- as.numeric(simulate(baseline, nsim = 10, innov = rep(0, 10)))
  expect_equal(x - y, rep(2, 10))
  for (seasonal in c(FALSE, TRUE)) {
    model <- arima_model(frequency = 4, p = 1, d = as.integer(!seasonal),
                         q = 0, P = 1, D = as.integer(seasonal), Q = 0,
                         phi = 0.5, Phi = 0.2, constant = 2, sigma = 1e-8)
    baseline <- arima_model(frequency = 4, p = 1, d = as.integer(!seasonal),
                            q = 0, P = 1, D = as.integer(seasonal), Q = 0,
                            phi = 0.5, Phi = 0.2, constant = 0, sigma = 1e-8)
    set.seed(1)
    x <- as.numeric(simulate(model, nsim = 10, innov = rep(0, 10)))
    set.seed(1)
    y <- as.numeric(simulate(baseline, nsim = 10, innov = rep(0, 10)))
    expect_equal(x - y, 2 * seq_len(10))
  }
  expect_error(arima_model(p = 0, d = 2, q = 0, constant = 1), "must be zero")
})

test_that("ETS rejects invalid fixed smoothing parameters promptly", {
  expect_error(ets_model(error = "A", trend = "N", seasonal = "N", alpha = 3),
               "between 0 and 1")
  expect_error(ets_model(error = "A", trend = "A", seasonal = "N",
                         damped = TRUE, phi = 0), "positive")
  expect_error(ets_model(error = "A", trend = "A", seasonal = "N",
                         alpha = 0.2, beta = 0.3), "beta must be <= alpha")
  expect_false(gratis:::ets_admissible(0.2, 0.1, NULL, 0, 1))
})

test_that("ETS resamples rejected draws and bounds retries", {
  draws <- numeric(0)
  fake_check <- function(alpha, ...) {
    draws <<- c(draws, alpha)
    length(draws) > 1
  }
  local_mocked_bindings(ets_admissible = fake_check, .package = "gratis")
  set.seed(1)
  model <- ets_model(error = "A", trend = "N", seasonal = "N")
  expect_length(draws, 2)
  expect_false(draws[1] == draws[2])
  expect_equal(unname(model$par["alpha"]), draws[2])
  local_mocked_bindings(ets_admissible = function(...) FALSE, .package = "gratis")
  expect_error(ets_model(error = "A", trend = "N", seasonal = "N", alpha = 0.2),
               "inadmissible")
  expect_error(ets_model(error = "A", trend = "N", seasonal = "N"),
               "could not find admissible")
})

test_that("multiplicative ETS initialization is positive and normalized", {
  for (m in c(4, 12, 365)) {
    set.seed(m)
    model <- ets_model(frequency = m, error = "M", trend = "N", seasonal = "M",
                       alpha = 0.2, gamma = 0.1)
    factors <- model$states[1, -1]
    expect_true(all(factors > 0))
    expect_equal(mean(factors), 1)
  }
  expect_error(ets_model(frequency = 4, seasonal = "M", season = c(1, 1, 1, -1)),
               "factors must be positive")
  expect_error(ets_model(error = "M", level = 0), "positive initial level")
})

test_that("GARCH standard deviations are converted to MAR variances", {
  skip_if_not_installed("fGarch")
  spec <- fGarch::garchSpec(model = list(omega = 16, alpha = 0, beta = 0))
  model <- mar_model(k = 1, p = 0, d = 0, constants = 0,
                     sigmas = list(spec), weights = 1)
  x <- simulate(model, nsim = 2000, seed = 1)
  expect_equal(sd(x), 4, tolerance = 0.2)
})
