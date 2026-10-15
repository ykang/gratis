test_that("repeated poles retain stationary variance after automatic burn-in", {
  phi <- c(5.4, -12.15, 14.58, -9.8415, 3.54294, -.531441)
  for (seasonal in c(FALSE, TRUE)) {
    model <- if (seasonal) {
      arima_model(frequency = 4, p = 0, d = 0, q = 0, P = 6, D = 0, Q = 0,
                  Phi = phi, constant = 0, sigma = 1)
    } else {
      arima_model(p = 6, d = 0, q = 0, phi = phi, constant = 0, sigma = 1)
    }
    burn <- gratis:::arima_burnin(model)
    impulse <- as.numeric(stats::filter(c(1, rep(0, 8191)), model$model$phi,
                                        method = "recursive"))
    ratio <- sum(impulse[seq_len(burn + 1)]^2) / sum(impulse^2)
    expect_gt(ratio, 1 - 1e-8)
    expect_gte(burn, 128)
  }
})

test_that("explicit ARIMA burn-in skips the energy search, not stationarity checks", {
  model <- arima_model(p = 1, d = 0, q = 0, phi = .5, constant = 0, sigma = 1)
  testthat::local_mocked_bindings(arima_burnin = function(...) stop("unexpected search"),
                                 .package = "gratis")
  expect_equal(as.numeric(simulate(model, innov = c(1, 0), n.start = 0)), c(1, .5))
  expect_equal(as.numeric(simulate(model, innov = 1, start.innov = 2)), 2)
  model$coef["ar1"] <- 1.1
  expect_error(simulate(model, nsim = 2, n.start = 0), "stationary")
})

test_that("unconverged ARIMA initialization fails with an explicit remedy", {
  model <- arima_model(p = 1, d = 0, q = 0, phi = .99999999, constant = 0, sigma = 1)
  expect_error(gratis:::arima_burnin(model), "supply n.start explicitly")
  expect_true(all(is.finite(simulate(model, nsim = 2, n.start = 10, seed = 1))))
})

test_that("R batch generation calculates automatic burn-in once", {
  model <- arima_model(p = 1, d = 0, q = 0, phi = .5, constant = 0, sigma = 1)
  original <- gratis:::arima_burnin
  calls <- 0L
  testthat::local_mocked_bindings(arima_burnin = function(object) {
    calls <<- calls + 1L
    original(object)
  }, .package = "gratis")
  result <- generate(model, length = 5, nseries = 3, seed = 7)
  expect_equal(calls, 1L)
  expect_equal(nrow(result), 15L)
  generate(model, length = 5, nseries = 3, seed = 7, n.start = 0)
  generate(model, length = 5, nseries = 3, seed = 7, start.innov = 1)
  expect_equal(calls, 1L)
})

test_that("near AR-MA cancellation does not conceal slow impulse energy", {
  phi <- .9999
  residual <- 1e-5
  model <- arima_model(p = 1, d = 0, q = 1, phi = phi, theta = residual - phi,
                        constant = 0, sigma = 1)
  burn <- gratis:::arima_burnin(model)
  variance <- 1 + residual^2 / (1 - phi^2)
  retained <- 1 + residual^2 * (1 - phi^(2 * burn)) / (1 - phi^2)
  expect_gt(retained / variance, 1 - 1e-8)
})

test_that("stationary AR sampling covers every possible real-root count", {
  set.seed(7)
  for (p in 3:6) {
    samples <- replicate(100, {
      phi <- gratis:::stationary_ar(p)
      roots <- polyroot(c(1, -phi))
      c(stationary = all(Mod(roots) > 1), count = sum(abs(Im(roots)) < 1e-7))
    })
    expect_true(all(samples["stationary", ] == 1))
    expect_equal(sort(unique(samples["count", ])), seq(p %% 2, p, by = 2))
  }
})

test_that("ARIMA and ETS reject unrepresentable noise variances", {
  for (sigma in c(1e155, 1e-200)) {
    expect_error(arima_model(sigma = sigma), "normal variance")
    expect_error(ets_model(sigma = sigma), "normal variance")
  }
  for (kind in c("arima", "ets")) {
    make_model <- if (kind == "arima") {
      function(sigma) arima_model(p = 0, d = 0, q = 0, constant = 0, sigma = sigma)
    } else {
      function(sigma) ets_model(error = "A", trend = "N", seasonal = "N",
                                alpha = .2, level = 0, sigma = sigma)
    }
    unit <- as.numeric(simulate(make_model(1), nsim = 20, seed = 7))
    for (sigma in c(1e150, 1e-150)) {
      scaled <- as.numeric(simulate(make_model(sigma), nsim = 20, seed = 7)) / sigma
      expect_equal(scaled, unit, tolerance = 1e-12)
    }
  }
})

test_that("fast MAR paths preserve recursion, short lengths and RNG draws", {
  cases <- list(c(0, 1), c(0, 3), c(1, 1), c(8, 1))
  for (case in cases) {
    p <- case[1]
    k <- case[2]
    phi <- matrix(0, p, k)
    if (p > 0) phi[p, ] <- .4
    model <- mar_model(k = k, phi = phi, d = 0, constants = seq_len(k) - .8,
                       sigmas = seq_len(k), weights = rep(1, k))
    for (size in list(c(1, 0), c(5, 0), c(30, 10))) {
      n <- size[1]
      burn <- size[2]
      total <- n + burn + p
      set.seed(7)
      components <- sample(seq_len(k), total, replace = TRUE, prob = model$weights)
      noise <- rnorm(total, sd = sqrt(model$sigmas[components]))
      expected <- numeric(total)
      for (i in seq.int(p + 1L, total)) {
        j <- components[i]
        expected[i] <- model$ar[1, j] + sum(model$ar[-1, j] * expected[i - seq_len(p)]) + noise[i]
      }
      expected_seed <- .Random.seed
      set.seed(7)
      actual <- simulate(model, nsim = n, n.start = burn)
      expect_equal(as.numeric(actual), tail(expected, n), tolerance = 1e-12)
      expect_identical(.Random.seed, expected_seed)
    }
  }
})

test_that("GA truncates oversized suggestions without changing their order", {
  suggestions <- matrix(seq(-1, 1, length.out = 40), ncol = 2)
  result <- gratis:::ga_ts(type = "real-valued", n = 1, min = c(-1, -1), max = c(1, 1),
                           fitness = function(par) list(value = -sum(par^2), x = par[1]),
                           suggestions = suggestions, popSize = 10, maxiter = 1,
                           monitor = FALSE, seed = 7)
  expect_equal(result@population, suggestions[1:10, ])
  expect_equal(as.numeric(result@popTs), suggestions[1:10, 1])
})

test_that("parallel updatePop keeps parameters, scores and paths synchronized", {
  skip_if_not_installed("doParallel")
  cluster <- try(parallel::makePSOCKcluster(2), silent = TRUE)
  if (inherits(cluster, "try-error")) skip("parallel worker sockets unavailable")
  on.exit({ parallel::stopCluster(cluster); foreach::registerDoSEQ() }, add = TRUE)
  fitness <- function(par) {
    value <- -.25^2
    attr(value, "par") <- .25
    list(value = value, x = .25)
  }
  arguments <- list(type = "real-valued", n = 1, min = -1, max = 1, fitness = fitness,
                    popSize = 10, maxiter = 2, run = 2, pcrossover = 0, pmutation = 0,
                    elitism = 0, updatePop = TRUE, monitor = FALSE, seed = 7)
  serial <- do.call(gratis:::ga_ts, arguments)
  result <- do.call(gratis:::ga_ts, c(arguments, list(parallel = cluster)))
  expect_equal(result@population, serial@population)
  expect_equal(result@population, matrix(.25, 10, 1))
  expect_equal(result@popTs, matrix(.25, 1, 10))
  expect_equal(result@fitness, rep(-.25^2, 10))
  expect_length(parallel::clusterCall(cluster, function() TRUE), 2)
})

test_that("target selection uses cached final fitness without storing history", {
  fake_ga <- function(...) {
    expect_false(list(...)$keepBest)
    methods::new("ga", popTs = cbind(c(1, 2, 3), c(4, 5, 6)), fitness = c(-.5, -.1))
  }
  testthat::local_mocked_bindings(ga_ts = fake_ga, .package = "gratis")
  result <- simulate_target(length = 3, target = 0, tolerance = .2,
                             feature_function = function(x) stop("must not recompute"))
  expect_equal(as.numeric(result), 4:6)
})

test_that("local GA optimization honors fitness-driven parameter updates", {
  for (poptim in c(0, 1)) {
    result <- gratis:::ga_ts(type = "real-valued", n = 1, min = -1, max = 1,
                             fitness = function(par) {
                               value <- -.25^2
                               attr(value, "par") <- .25
                               list(value = value, x = .25)
                             }, updatePop = TRUE, optim = TRUE,
                             optimArgs = list(poptim = poptim),
                             popSize = 10, maxiter = 1, monitor = FALSE, seed = 7)
    expect_equal(result@population, matrix(.25, 10, 1))
    expect_equal(result@popTs, matrix(.25, 1, 10))
    expect_equal(result@fitness, rep(-.25^2, 10))
  }
})

test_that("Monte Carlo MAR moments match their analytic DGP", {
  for (mixture in c(FALSE, TRUE)) {
    a <- if (mixture) c(.4, -.2) else .6
    constant <- if (mixture) c(.5, -.2) else .2
    v <- if (mixture) c(1, 2) else 1
    w <- if (mixture) c(.7, .3) else 1
    mu <- sum(w * constant) / (1 - sum(w * a))
    second <- (sum(w * (constant^2 + v)) + 2 * sum(w * constant * a) * mu) /
      (1 - sum(w * a^2))
    expected <- c(mu, second - mu^2, sum(w * a))
    model <- mar_model(phi = matrix(a, nrow = 1), d = 0, constants = constant,
                       sigmas = v, weights = w)
    estimates <- vapply(seq_len(20), function(seed) {
      x <- as.numeric(simulate(model, nsim = 2000, n.start = 300, seed = seed))
      c(mean(x), var(x), cor(head(x, -1), tail(x, -1)))
    }, numeric(3))
    mcse <- apply(estimates, 1, sd) / sqrt(NCOL(estimates))
    expect_true(all(abs(rowMeans(estimates) - expected) < 6 * mcse + .005))
  }
})
