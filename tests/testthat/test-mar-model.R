test_that("mar_model handles zero-order zero-constant models", {
  model <- mar_model(
    k = 1,
    p = 0,
    d = 0,
    constants = 0,
    sigmas = 1,
    weights = 1
  )

  expect_s3_class(model, "mar")
  expect_equal(dim(model$ar), c(1L, 1L))
  expect_equal(unname(model$ar[1, 1]), 0)
})

test_that("mar_model uses supplied seasonal coefficient matrix order", {
  phi <- matrix(0.1, nrow = 1, ncol = 2)
  Phi <- matrix(c(0.2, 0.1, 0.3, 0.2), nrow = 2, ncol = 2)

  model <- mar_model(
    phi = phi,
    Phi = Phi,
    d = 0,
    D = 0,
    sigmas = c(1, 1),
    weights = c(0.5, 0.5),
    seasonal_periods = 12
  )

  expect_equal(model$P, c(2L, 2L))
  expect_equal(model$Phi, Phi)
})

test_that("mar_model validates seasonal period dimensions", {
  expect_error(
    mar_model(k = 2, seasonal_periods = c(12, 24, 168)),
    "Dimension of seasonal_periods"
  )
})
