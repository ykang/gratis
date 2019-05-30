generate_ts_with_target_ts <- function(n, ts.length, freq, seasonal, x0, max.fitness = -3, h = 8, parallel=TRUE) {
  ga_min <-
    if (seasonal == 0) {
      c(rep(-0.5, 6), rep(0, 4))
    } else if (seasonal == 1) {
      c(rep(ifelse(nsdiffs(x0) == 1, -0.5, -1), 12), rep(0, 5))
    } else {
      c(rep(0, 35))
    }
  ga_max <-
    if (seasonal == 0) {
      c(rep(0.5, 6), rep(1, 4))
    } else if (seasonal == 1) {
      c(rep(ifelse(nsdiffs(x0) == 1, 0.5, 1), 12), rep(0.5, 5))
    } else {
      c(rep(1, 35))
    }
  evolved.ts <- c()
  while (ifelse(is.null(dim(evolved.ts)), 0 < 1, dim(evolved.ts)[2] < n)) {
    GA <- ga_ts(
      type = "real-valued", fitness = fitness_ts1, x0 = x0, seasonal = seasonal,
      ts.length, freq, 3, h = h,
      n = ts.length,
      min = ga_min,
      max = ga_max,
      parallel = parallel, popSize = 30, maxiter = 100,
      pmutation = 0.3, pcrossover = 0.8, maxFitness = max.fitness,
      run = 10, keepBest = TRUE, monitor = GA::gaMonitor
    )
    evolved.ts.new <-
      unique(do.call(
        cbind,
        eval(parse(text = paste("list(", paste("GA@bestSol[[GA@iter - ", 0:(GA@run - 1), "]]", sep = "", collapse = ","), ")")))
      ), MARGIN = 2)
    evolved.ts <- cbind(evolved.ts, evolved.ts.new)
  }
  if (length(freq) == 1) {
    evolved.ts <- ts(evolved.ts[, 1:n], frequency = freq)
  } else {
    evolved.ts <- msts(evolved.ts[, 1:n], seasonal.periods = freq)
  }
  return(evolved.ts)
}
