test_that("simulate_target honors user-supplied component count", {
  captured <- NULL
  fake_ga <- function(...) {
    captured <<- list(...)
    stop("stop after capturing ga_ts arguments", call. = FALSE)
  }

  testthat::local_mocked_bindings(ga_ts = fake_ga, .package = "gratis")

  expect_error(
    simulate_target(
      length = 20,
      feature_function = function(x) mean(x),
      target = 0,
      k = 2
    ),
    "stop after capturing ga_ts arguments",
    fixed = TRUE
  )

  expect_equal(captured$ncomponents, 2)
  expect_length(captured$min, 15)
  expect_length(captured$max, 15)
})
