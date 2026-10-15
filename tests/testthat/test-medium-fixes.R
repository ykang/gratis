test_that("Shiny is no longer a package API or suggested dependency", {
  expect_false("app_gratis" %in% getNamespaceExports("gratis"))
  expect_false(exists("app_gratis", envir = asNamespace("gratis"), inherits = FALSE))
  suggested <- strsplit(packageDescription("gratis")$Suggests, ",\\s*")[[1]]
  expect_false(any(c("shiny", "shinydashboard", "rlang") %in% trimws(suggested)))
})

batch_models <- function() {
  list(mar_model(k = 1, p = 1, phi = .4, d = 0, constants = 0,
                 sigmas = 1, weights = 1, seasonal_periods = 4),
       arima_model(frequency = 4, p = 1, d = 0, q = 0, P = 0, D = 0, Q = 0,
                   phi = .4, constant = 0, sigma = 1),
       ets_model(frequency = 4, error = "A", trend = "N", seasonal = "N",
                 alpha = .2, level = 0, sigma = 1))
}

test_that("R batch seeds generate distinct reproducible paths and restore the RNG", {
  for (model in batch_models()) {
    set.seed(123)
    previous <- .Random.seed
    first <- generate(model, length = 12, nseries = 2, seed = 7)
    expect_identical(.Random.seed, previous)
    expect_identical(first, generate(model, length = 12, nseries = 2, seed = 7))
    paths <- split(first$value, first$key)
    expect_false(identical(paths[[1]], paths[[2]]))
    set.seed(7)
    reference <- generate(model, length = 12, nseries = 2)
    expect_identical(first, reference)
    expect_false(identical(.Random.seed, previous))
  }
})

test_that("batch seeding restores an absent RNG state and restores state on errors", {
  model <- batch_models()[[1]]
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous <- .Random.seed
  on.exit(if (had_seed) {
    assign(".Random.seed", previous, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(generate(model, length = 12, nseries = 2, seed = 7))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  set.seed(123)
  previous_error <- .Random.seed
  expect_error(generate(model, length = 12, nseries = 2, seed = 7, n.start = -1), "non-negative")
  expect_identical(.Random.seed, previous_error)
})

test_that("target generation preserves the same time index as model generation", {
  local_mocked_bindings(simulate_target = function(length, seasonal_periods, ...) {
    ts(seq_len(length), frequency = max(seasonal_periods))
  }, .package = "gratis")
  for (period in c(1, 4, 12, 24, 52)) {
    model <- mar_model(k = 1, p = 0, d = 0, constants = 0, sigmas = 1,
                       weights = 1, seasonal_periods = period)
    reference <- generate(model, length = 12, nseries = 2, seed = 7)
    actual <- generate_target(length = 12, nseries = 2, seasonal_periods = period,
                              feature_function = mean, target = 0)
    expect_identical(actual$index, reference$index)
    expect_identical(actual$key, reference$key)
  }
})

test_that("nonuniform mutation moves in both directions within bounds", {
  object <- gratis:::ga_ts(type = "real-valued", n = 1,
                           fitness = function(p) list(value = -p[1]^2, x = p[1]),
                           min = -2, max = 4, suggestions = matrix(1, 10, 1),
                           popSize = 10, maxiter = 1, monitor = FALSE)
  object@maxiter <- 100L
  set.seed(1)
  mutations <- replicate(200, gratis:::gareal_nraMutation(object, 1))
  expect_true(any(mutations < 1))
  expect_true(any(mutations > 1))
  expect_true(all(mutations >= -2 & mutations <= 4))
  object@iter <- object@maxiter
  expect_equal(gratis:::gareal_nraMutation(object, 1), 1)
})

test_that("permutation initialization respects arbitrary integer bounds", {
  for (bounds in list(c(3, 6), c(-2, 2), c(3, 3))) {
    expected <- seq(bounds[1], bounds[2])
    object <- gratis:::ga_ts(type = "permutation", n = 1,
                             fitness = function(p) list(value = -sum(p), x = p[1]),
                             min = bounds[1], max = bounds[2], popSize = 10,
                             maxiter = 1, monitor = FALSE, seed = 1)
    expect_equal(dim(object@population), c(10L, length(expected)))
    expect_true(all(apply(object@population, 1, function(p) identical(sort(p), as.numeric(expected)))))
  }
})

test_that("vectorized conditional densities match pointwise references", {
  parameters <- list(c(.1, .2), c(-.1, .4, .1), c(0, .6))
  weights <- c(.2, .3, .5)
  for (n in c(0, 1, 2, 3, 100)) {
    set.seed(1)
    y <- if (n) ts(rnorm(n)) else numeric(0)
    variances <- list(seq(.5, 2, length.out = n), rep(1, n), rep(2, n))
    expected <- numeric(n)
    for (i in seq_len(max(0, n - 2)) + 2L) {
      means <- vapply(parameters, function(p) sum(p * c(1, y[i - seq_len(length(p) - 1L)])), numeric(1))
      expected[i] <- gratis:::dmixnorm(y[i], means, vapply(variances, `[`, numeric(1), i), weights, log = TRUE)
    }
    actual <- gratis:::dmixnorm_ts(y, parameters, variances, weights, log = TRUE)
    expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
    if (n) expect_s3_class(actual, "ts") else expect_length(actual, 0)
    expect_equal(as.numeric(gratis:::dmixnorm_ts(y, parameters, variances, weights)), exp(expected))
  }
})

test_that("inactive components cannot poison R conditional densities", {
  actual <- gratis:::dmixnorm_ts(c(1e308, 1e308, 0), list(c(0, 0, 0), c(0, 2, -2)),
                                 list(rep(1, 3), rep(1, 3)), c(1, 0), log = TRUE)
  expect_equal(actual, c(0, 0, dnorm(0, log = TRUE)))
})
