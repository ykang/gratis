test_that("MAR inputs reject nonfinite values and fractional orders", {
  base <- list(k = 1, p = 0, d = 0, constants = 0, sigmas = 1, weights = 1)
  for (name in c("phi", "Phi", "constants", "sigmas")) {
    for (value in c(NA_real_, Inf, -Inf)) {
      arguments <- base
      arguments[[name]] <- value
      expect_error(do.call(mar_model, arguments), "finite")
    }
  }
  for (name in c("p", "d", "P", "D", "k", "seasonal_periods")) {
    for (value in c(1.9, NA_real_, Inf, -1)) {
      arguments <- base
      arguments[[name]] <- value
      expect_error(do.call(mar_model, arguments), "integer")
    }
  }
  for (name in c("ar", "ma", "sar", "sma")) {
    expect_error(do.call(pi_coefficients, setNames(list(NA_real_), name)), "finite")
  }
})

test_that("R mixture weights and tail log densities follow probability conventions", {
  expect_equal(as.numeric(gratis:::dmixnorm(0, 0, 1, 2)), dnorm(0))
  expect_equal(as.numeric(gratis:::dmixnorm(40, 0, 1, 1, log = TRUE)), dnorm(40, log = TRUE))
  expect_equal(as.numeric(gratis:::dmixnorm(c(0, 40), c(0, 1), c(1, 2), c(1e308, 1e308), log = TRUE)),
               vapply(c(0, 40), function(x) {
                 values <- c(dnorm(x, 0, 1, log = TRUE), dnorm(x, 1, sqrt(2), log = TRUE))
                 max(values) + log(sum(exp(values - max(values)))) - log(2)
               }, numeric(1)))
  expect_identical(as.numeric(gratis:::dmixnorm(Inf, 0, 1, 1, log = TRUE)), -Inf)
  expect_length(gratis:::dmixnorm(numeric(), 0, 1, 1), 0)
})

test_that("R mixture sampling and density reject invalid covariance and weights", {
  invalid <- list(-diag(2), matrix(c(1, 2, 2, 1), 2), diag(c(1, 0)),
                  matrix(c(1, 0, 1, 1), 2), diag(c(NA_real_, 1)), diag(c(Inf, 1)))
  for (covariance in invalid) {
    arguments <- list(means = matrix(0, 2, 1), sigmas = array(covariance, c(2, 2, 1)), weights = 1)
    expect_error(do.call(rmixnorm, c(list(n = 2), arguments)), "positive definite")
    expect_error(do.call(gratis:::dmixnorm, c(list(x = c(0, 0)), arguments)), "positive definite")
  }
  for (weights in list(c(0, 0), c(-1, 1), c(NA_real_, 1), c(Inf, 1))) {
    expect_error(rmixnorm(1, c(0, 1), c(1, 1), weights), "weights")
    expect_error(gratis:::dmixnorm(0, c(0, 1), c(1, 1), weights), "weights")
  }
  expect_error(rmixnorm(0, c(0, 0), c(1, -1), c(1, 0)), "positive definite")
})

test_that("sparse MAR recursion preserves the dense seeded path", {
  for (period in c(1, 12, 365)) {
    model <- mar_model(phi = matrix(c(.2, .1), 1), d = 0, P = 2, D = 0,
                       Phi = matrix(c(.3, .1, .5, -.1), 2), constants = c(.1, -.1),
                       sigmas = c(1, 2), weights = c(.4, .6), seasonal_periods = period)
    set.seed(7)
    p <- NROW(model$ar) - 1L
    n <- 300 + 100 + p
    components <- sample(seq_len(2), n, replace = TRUE, prob = model$weights)
    noise <- rnorm(n, sd = sqrt(model$sigmas[components]))
    expected <- numeric(n)
    for (i in (p + 1L):n) {
      expected[i] <- sum(model$ar[, components[i]] * c(1, expected[i - seq_len(p)])) + noise[i]
    }
    expect_equal(as.numeric(simulate(model, 300, seed = 7)), tail(expected, 300), tolerance = 1e-12)
  }
})

test_that("scalar mixture recursion preserves general sampling and density", {
  parameters <- list(c(.1, .2), c(-.1, .4, .1), c(0, .6))
  variances <- list(seq(.5, 2, length.out = 100), rep(1, 100), rep(2, 100))
  weights <- c(.2, .3, .5)
  expected <- numeric(100)
  expected_log <- numeric(100)
  set.seed(7)
  for (i in 3:100) {
    means <- vapply(parameters, function(p) sum(p * c(1, expected[i - seq_len(length(p) - 1L)])), numeric(1))
    sigmas <- vapply(variances, `[`, numeric(1), i)
    expected[i] <- rmixnorm(1, means, sigmas, weights)
    expected_log[i] <- gratis:::dmixnorm(expected[i], means, sigmas, weights, log = TRUE)
  }
  set.seed(7)
  actual <- rmixnorm_ts(100, parameters, variances, weights)
  expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
  expect_equal(as.numeric(gratis:::dmixnorm_ts(actual, parameters, variances, weights, log = TRUE)),
               expected_log, tolerance = 1e-12)
  for (variance in c(0, -1, NA_real_, Inf)) {
    expect_error(rmixnorm_ts(1, list(c(0, .2), c(0, .5)), list(1, variance), c(1, 0)), "positive and finite")
    expect_error(gratis:::dmixnorm_ts(0, list(c(0, .2), c(0, .5)), list(1, variance), c(1, 0)),
                 "positive and finite")
  }
})

test_that("parallel GA preserves scalar series and caller-owned clusters", {
  skip_if_not_installed("doParallel")
  cluster <- try(parallel::makePSOCKcluster(2), silent = TRUE)
  if (inherits(cluster, "try-error")) skip("parallel worker sockets unavailable")
  on.exit({ parallel::stopCluster(cluster); foreach::registerDoSEQ() }, add = TRUE)
  result <- gratis:::ga_ts(type = "real-valued", n = 1, min = -1, max = 1,
                           fitness = function(par) list(value = -par[1]^2, x = par[1]),
                           parallel = cluster, popSize = 10, maxiter = 2, run = 2,
                           pcrossover = 0, pmutation = 0, elitism = 0, monitor = FALSE, seed = 7)
  expect_equal(dim(result@popTs), c(1L, 10L))
  expect_equal(unname(result@popTs[1, ]), unname(result@population[, 1]))
  expect_equal(result@fitness, -result@population[, 1]^2, ignore_attr = TRUE)
  expect_length(parallel::clusterCall(cluster, function() TRUE), 2L)
})
