test_that("rmixnorm supports any number of mixture components", {
  set.seed(1)

  x <- rmixnorm(
    n = 200,
    means = seq_len(6),
    sigmas = rep(1, 6),
    weights = rep(1 / 6, 6)
  )

  expect_type(x, "double")
  expect_length(x, 200)
  expect_false(anyNA(x))
})

test_that("dmixnorm_ts returns densities consistent with log densities", {
  y <- ts(c(0, 0.1, -0.2, 0.3))
  means <- list(c(0, 0.5), c(0, -0.2))
  sigmas <- list(rep(1, length(y)), rep(2, length(y)))
  weights <- c(0.4, 0.6)

  log_density <- gratis:::dmixnorm_ts(
    y = y,
    means.ar.par.list = means,
    sigmas.list = sigmas,
    weights = weights,
    log = TRUE
  )
  density <- gratis:::dmixnorm_ts(
    y = y,
    means.ar.par.list = means,
    sigmas.list = sigmas,
    weights = weights,
    log = FALSE
  )

  expect_equal(as.numeric(density), exp(as.numeric(log_density)))
  expect_true(all(as.numeric(density) > 0))
})
