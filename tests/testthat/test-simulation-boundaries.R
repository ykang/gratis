test_that("MAR simulation preserves lengths at zero burn-in and zero order", {
  model <- mar_model(k = 1, p = 0, d = 0, constants = 0, sigmas = 1, weights = 1)
  for (n in c(1, 5)) expect_length(simulate(model, nsim = n, n.start = 0), n)
  expect_length(simulate(model, nsim = 0, n.start = 0), 0)
  expect_error(simulate(model, nsim = -1), "non-negative integers")
  expect_error(simulate(model, n.start = 1.5), "non-negative integers")
})

test_that("short mixture series return their initialization without descending loops", {
  for (n in 0:3) {
    pars <- list(c(0, .5, .2, .1))
    variances <- list(rep(1, n))
    x <- rmixnorm_ts(n, pars, variances, 1, yinit = 2)
    expect_equal(as.numeric(x), rep(2, n))
    expect_equal(as.numeric(gratis:::dmixnorm_ts(x, pars, variances, 1)), rep(1, n))
    expect_equal(as.numeric(gratis:::dmixnorm_ts(x, pars, variances, 1, log = TRUE)), rep(0, n))
  }
})

test_that("MAR supports zero and large finite weights but rejects invalid weights", {
  args <- list(k = 2, p = 0, d = 0, constants = c(0, 0), sigmas = c(1, 1))
  expect_equal(do.call(mar_model, c(args, list(weights = c(1, 0))))$weights, c(1, 0))
  expect_equal(do.call(mar_model, c(args, list(weights = c(1e308, 1e308))))$weights, c(.5, .5))
  for (weights in list(c(0, 0), c(NA, 1), c(Inf, 1), c(-1, 1))) {
    expect_error(do.call(mar_model, c(args, list(weights = weights))), "weights must")
  }
})
