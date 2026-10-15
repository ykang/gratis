test_that("ARIMA simulation scales output and burn-in innovations together", {
  unit <- arima_model(p = 1, d = 0, q = 1, phi = 0.99, theta = 0.3,
                      constant = 0, sigma = 1)
  tiny <- unit
  tiny$sigma2 <- 1e-16
  expect_s3_class(unit, "gratis_arima")
  expect_s3_class(unit, "Arima")
  x <- as.numeric(simulate(unit, nsim = 100, seed = 1))
  y <- as.numeric(simulate(tiny, nsim = 100, seed = 1))
  expect_equal(y, x * 1e-8, tolerance = 1e-12)
  expect_lt(max(abs(y)), 1e-5)
  set.seed(123)
  previous_seed <- .Random.seed
  simulate(unit, nsim = 5, seed = 1)
  expect_identical(.Random.seed, previous_seed)
})

test_that("ARIMA automatic burn-in accounts for near-unit roots", {
  phi <- 0.999
  sigma <- 0.01
  model <- arima_model(p = 1, d = 0, q = 0, phi = phi, constant = 0, sigma = sigma)
  burn <- gratis:::arima_burnin(model)
  expect_lt(phi^(2 * (burn + 1L)), 1e-8)
  set.seed(1)
  noise <- rnorm(burn + 3L, sd = sigma)
  expected <- vapply(burn + seq_len(3L), function(t) {
    sum(noise[seq_len(t)] * phi^(t - seq_len(t)))
  }, numeric(1))
  expect_equal(as.numeric(simulate(model, nsim = 3, seed = 1)), expected,
               tolerance = 1e-10)
})

test_that("ARIMA preserves zero-sum seasonal AR and MA coefficients", {
  ar <- arima_model(frequency = 4, p = 0, d = 0, q = 0, P = 2, D = 0, Q = 0,
                    Phi = c(0.5, -0.5), constant = 0, sigma = 1)
  ma <- arima_model(frequency = 4, p = 0, d = 0, q = 0, P = 0, D = 0, Q = 2,
                    Theta = c(0.5, -0.5), constant = 0, sigma = 1)
  innovations <- c(1, rep(0, 11))
  expected_ar <- expected_ma <- numeric(12)
  expected_ar[c(1, 5, 9)] <- c(1, 0.5, -0.25)
  expected_ma[c(1, 5, 9)] <- c(1, 0.5, -0.5)
  expect_equal(as.numeric(simulate(ar, innov = innovations, n.start = 0)), expected_ar)
  expect_equal(as.numeric(simulate(ma, innov = innovations, n.start = 0)), expected_ma)
  expect_equal(as.numeric(simulate(ar, innov = 1, n.start = 0)), 1)
  expect_equal(as.numeric(simulate(ma, innov = 1, n.start = 0)), 1)
})

test_that("ARIMA supports controlled innovations and explicit zero states", {
  model <- arima_model(p = 1, d = 0, q = 1, phi = 0.5, theta = 0.2,
                       constant = 0, sigma = 2)
  innovations <- c(0.1, -0.3, 0.2, 0.5)
  expect_equal(as.numeric(simulate(model, innov = innovations, n.start = 0)),
               c(0.1, -0.23, 0.025, 0.5525))
  expect_equal(as.numeric(simulate(model, innov = innovations, start.innov = 0.4)),
               c(0.38, -0.09, 0.095, 0.5875))
  expect_error(simulate(model, n.start = -1), "non-negative integer")
  expect_error(simulate(model, n.start = 1.5), "non-negative integer")
  expect_error(simulate(model, innov = NA_real_), "finite numeric vector")
  expect_error(simulate(model, start.innov = Inf), "finite numeric vector")
  expect_error(simulate(model, n.start = 2, start.innov = 1), "length n.start")
})

test_that("generate dispatches package ARIMA models without changing fitted models", {
  model <- arima_model(p = 1, d = 0, q = 0, phi = 0.5, constant = 0, sigma = 1)
  generated <- generate(model, length = 12, nseries = 2)
  expect_s3_class(generated, "tbl_ts")
  expect_equal(nrow(generated), 24L)
  fitted <- forecast::Arima(ts(seq_len(30)), order = c(0, 0, 0))
  expect_false(inherits(fitted, "gratis_arima"))
})
