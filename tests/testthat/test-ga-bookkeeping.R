ga_test_fitness <- function(par, center = .25) {
  list(value = -(par[1] - center)^2, x = rep(par[1], 3))
}

check_ga_bookkeeping <- function(result) {
  finite <- is.finite(result@fitness)
  expect_equal(unname(result@popTs[, finite, drop = FALSE]),
               matrix(rep(result@population[finite, 1], each = 3), nrow = 3))
  expect_equal(unname(result@fitness[finite]), -(result@population[finite, 1] - .25)^2)
}

test_that("all real-valued selection operators preserve population time series", {
  selectors <- c("ga_lrSelection", "ga_nlrSelection", "ga_rwSelection", "ga_tourSelection",
                 "gareal_lsSelection", "gareal_sigmaSelection")
  for (selector in selectors) {
    result <- gratis:::ga_ts(type = "real-valued", fitness = ga_test_fitness, n = 3,
                            min = -1, max = 1, popSize = 10, maxiter = 2, run = 2,
                            pcrossover = 0, pmutation = 0, elitism = 0,
                            selection = getFromNamespace(selector, "gratis"),
                            keepBest = TRUE, monitor = FALSE, seed = 1)
    check_ga_bookkeeping(result)
    expect_true(all(is.finite(result@popTs)))
  }
})

test_that("proportional selectors handle invalid and equal fitness values", {
  for (selector in c("gareal_lsSelection", "gareal_sigmaSelection", "ga_rwSelection")) {
    result <- gratis:::ga_ts(type = "real-valued", n = 3, min = -1, max = 1,
                            fitness = function(par) {
                              if (par[1] < 0) list(value = -Inf, x = rep(NA_real_, 3))
                              else ga_test_fitness(par)
                            }, suggestions = matrix(c(rep(-1, 9), .5), ncol = 1),
                            popSize = 10, maxiter = 2, run = 2, pcrossover = 0,
                            pmutation = 0, elitism = 0, monitor = FALSE, seed = 1,
                            selection = getFromNamespace(selector, "gratis"))
    expect_true(all(is.finite(result@fitness)))
    check_ga_bookkeeping(result)
  }
})

test_that("local optimization updates scalar fitness, parameters and stored series", {
  for (poptim in c(0, 1)) {
    result <- gratis:::ga_ts(type = "real-valued", fitness = ga_test_fitness, center = .25,
                            n = 3, min = -1, max = 1, suggestions = matrix(.8, 10, 1),
                            popSize = 10, maxiter = 1, optim = TRUE,
                            optimArgs = list(poptim = poptim, control = list(maxit = 50)),
                            keepBest = TRUE, monitor = FALSE, seed = 1)
    expect_gt(result@fitnessValue, -1e-10)
    check_ga_bookkeeping(result)
    expect_true(any(abs(result@bestSol[[1]][1, ] - .25) < 1e-5))
  }
})

test_that("crossover invalidates cached fitness together with its series", {
  crossover <- function(object, parents) {
    list(children = object@population[rev(parents), , drop = FALSE],
         fitness = object@fitness[rev(parents)])
  }
  result <- gratis:::ga_ts(type = "real-valued", fitness = ga_test_fitness, n = 3,
                          min = -1, max = 1, popSize = 10, maxiter = 2, run = 2,
                          pcrossover = 1, crossover = crossover, pmutation = 0,
                          elitism = 0, keepBest = TRUE, monitor = FALSE, seed = 1)
  expect_true(all(is.finite(result@popTs)))
  check_ga_bookkeeping(result)
})

test_that("binary and permutation selection aliases preserve series bookkeeping", {
  for (type in c("binary", "permutation")) {
    args <- list(type = type, fitness = ga_test_fitness, n = 3, popSize = 10,
                 maxiter = 2, run = 2, keepBest = TRUE, monitor = FALSE, seed = 1)
    if (type == "binary") args$nBits <- 4
    else {
      args$min <- 1
      args$max <- 4
    }
    result <- do.call(gratis:::ga_ts, args)
    check_ga_bookkeeping(result)
    expect_true(all(is.finite(result@popTs)))
  }
})

test_that("elitism handles fewer unique individuals than requested elites", {
  result <- gratis:::ga_ts(type = "real-valued", fitness = ga_test_fitness, n = 3,
                          min = -1, max = 1, suggestions = matrix(.8, 10, 1),
                          popSize = 10, maxiter = 2, run = 2, pcrossover = 0,
                          pmutation = 0, elitism = 5, monitor = FALSE, seed = 1)
  check_ga_bookkeeping(result)
  expect_true(all(is.finite(result@population)))
})

test_that("post-fitness updates are reflected in the summary and saved best series", {
  post <- function(object, ...) {
    object@population[1, 1] <- .25
    object@fitness[1] <- 0
    object@popTs[, 1] <- .25
    object
  }
  result <- gratis:::ga_ts(type = "real-valued", fitness = ga_test_fitness, n = 3,
                          min = -1, max = 1, suggestions = matrix(.8, 10, 1),
                          popSize = 10, maxiter = 2, postFitness = post,
                          maxFitness = -1e-10, keepBest = TRUE, monitor = FALSE, seed = 1)
  expect_equal(result@iter, 1)
  expect_equal(unname(result@summary[1, "max"]), 0)
  expect_equal(unname(result@bestSol[[1]]), matrix(.25, 3, 1))
  check_ga_bookkeeping(result)
})
