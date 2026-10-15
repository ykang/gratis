target_parameters <- function(weights) {
  k <- length(weights)
  c(rep(0, 6 * k), rep(1, k), head(cumsum(weights), -1L))
}

test_that("target weights cover the full simplex, including one component", {
  original_mar <- mar_model
  captured <- NULL
  testthat::local_mocked_bindings(mar_model = function(...) {
    model <- original_mar(...)
    captured <<- model$weights
    model
  }, .package = "gratis")
  for (weights in list(1, c(.9, .05, .05), c(.05, .9, .05), c(0, 0, 1),
                      c(1, 0, 0), c(.2, 0, .8))) {
    fit <- gratis:::fitness_mar(target_parameters(weights), length = 12,
                                seasonal_periods = 1, ncomponents = length(weights),
                                feature_function = mean, target = 0)
    expect_equal(captured, weights)
    expect_true(is.finite(fit$value))
    expect_true(all(is.finite(fit$x)))
  }
})

test_that("scaling is stable and preserves time-series attributes", {
  inputs <- list(c(1e155, 2e155, 3e155), 1000001:1000003, c(-1e308, 0, 1e308),
                 c(1e-300, 2e-300, 3e-300), 1 + (0:2) * .Machine$double.eps)
  for (x in inputs) {
    expect_equal(gratis:::scalets(x), c(-1, 0, 1), tolerance = 1e-12)
  }
  x <- forecast::msts(1000001:1000012, seasonal.periods = c(4, 12), start = 2000)
  scaled <- gratis:::scalets(x)
  expect_identical(attributes(scaled), attributes(x))
  expect_equal(mean(scaled), 0, tolerance = 1e-12)
  expect_equal(sd(scaled), 1)
  for (x in list(numeric(0), 5, c(1e155, 1e155), c(0, 0))) {
    expect_equal(gratis:::scalets(x), x)
  }
  expect_error(gratis:::scalets(c(Inf, 1)), "finite")
  expect_error(gratis:::scalets(c(NA_real_, 1)), "finite")
})

test_that("explosive candidates cannot become perfect zero-series matches", {
  fit <- gratis:::fitness_mar(c(0, 0, 1.5, 0, 0, 0, 1), length = 1000,
                              seasonal_periods = 1, ncomponents = 1,
                              feature_function = sd, target = 0)
  expect_equal(fit$value, -1)
  expect_equal(sd(fit$x), 1)
  fit <- gratis:::fitness_mar(c(0, 0, 1.5, 0, 0, 0, 1), length = 4000,
                              seasonal_periods = 1, ncomponents = 1,
                              feature_function = function(x) stop("must not be called"), target = 0)
  expect_identical(fit$value, -Inf)
})

test_that("target search reports tolerance misses and callback failures", {
  original_ga <- gratis:::ga_ts
  fast_ga <- function(...) {
    args <- list(...)
    args$popSize <- 10
    args$maxiter <- 2
    args$run <- 2
    do.call(original_ga, args)
  }
  testthat::local_mocked_bindings(ga_ts = fast_ga, .package = "gratis")
  expect_warning(x <- simulate_target(length = 20, k = 1, feature_function = function(x) 0,
                                      target = 1, tolerance = .01),
                 "mean absolute feature error=1, tolerance=0.01")
  expect_length(x, 20)
  expect_true(all(is.finite(x)))
  x <- simulate_target(length = 20, k = 1, feature_function = mean, target = 0, tolerance = 1e-8)
  expect_lte(abs(mean(x)), 1e-8)
  expect_error(simulate_target(length = 20, k = 1,
                               feature_function = function(x) stop("callback bug"), target = 0),
               "callback bug")
  expect_error(simulate_target(length = 20, k = 1, feature_function = function(x) c(0, 1), target = 0),
               "one numeric value per target feature")
  expect_error(simulate_target(length = 20, k = 1, feature_function = function(x) NA_real_, target = 0),
               "finite value")
})

test_that("final target selection uses mean absolute error, not squared error", {
  fake_ga <- function(...) {
    args <- list(...)
    expect_false(args$keepBest)
    candidates <- cbind(c(.2, 0, 0, 0), c(.11, .11, 0, 0))
    scores <- apply(candidates, 2, gratis:::target_score,
                    feature_function = args$feature_function, target = args$target)
    methods::new("ga", popTs = candidates, fitness = scores)
  }
  testthat::local_mocked_bindings(ga_ts = fake_ga, .package = "gratis")
  x <- simulate_target(length = 4, feature_function = function(x) as.numeric(x[1:2]),
                       target = c(0, 0), tolerance = .105)
  expect_equal(as.numeric(x), c(.2, 0, 0, 0))
})

test_that("target search validates targets and tolerances", {
  expect_error(simulate_target(feature_function = mean, target = numeric(0)), "nonempty finite")
  expect_error(simulate_target(feature_function = mean, target = NA_real_), "nonempty finite")
  expect_error(simulate_target(feature_function = mean, target = 0, tolerance = -1), "tolerance")
  expect_error(simulate_target(feature_function = mean, target = 0, tolerance = Inf), "tolerance")
})
