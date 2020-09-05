#' @importFrom forecast nsdiffs
#' @importFrom forecast BoxCox
#' @importFrom forecast InvBoxCox
#' @importFrom tsibble as_tsibble

generate_ts_with_target_ts <- function(n, ts.length, freq, seasonal, x, max.fitness = -3, h = 8, preprocessing = 1, parallel=TRUE, output_format="list") {
  ga_min <-
    if (seasonal == 0) {
      c(rep(-0.5, 6), rep(0, 4))
    } else if (seasonal == 1) {
      c(rep(ifelse(nsdiffs(x) == 1, -0.5, -1), 12), rep(0, 5))
    } else {
      c(rep(0, 35))
    }
  ga_max <-
    if (seasonal == 0) {
      c(rep(0.5, 6), rep(1, 4))
    } else if (seasonal == 1) {
      c(rep(ifelse(nsdiffs(x) == 1, 0.5, 1), 12), rep(0.5, 5))
    } else {
      c(rep(1, 35))
    }

  if (preprocessing==1){
    out_req <- Smoothing_ts2(x, h, h)
    x0 <- out_req$series
    SeasIndIn <- out_req$seasonalIn
    SeasInd <- out_req$seasonal
    lambda <- out_req$lambda
  }else{
    x0 <- x
    SeasIndIn <- rep(0, length(x))
    SeasInd <- rep(0, h)
  }
  x0.mean = mean(x0)
  x0.sd = sd(x0)
  x0 <- (x0 - x0.mean)/x0.sd
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
    evolved.ts <- ts(x0.sd * evolved.ts[, 1:n] + x0.mean, frequency = freq, start = attributes(x)$tsp[1])
  } else {
    evolved.ts <- msts(evolved.ts[, 1:n], seasonal.periods = freq)
  }
  if (all(SeasIndIn==0)==F){
    for (i in 1:n){
      evolved.ts[,i] <- ts(InvBoxCox((BoxCox(evolved.ts[,i], lambda) + c(SeasIndIn, SeasInd)) , lambda),
                              frequency = freq, start = attributes(x)$tsp[1])
    }
  }
  
  # New content 
  output <- if (output_format == "list") {
    evolved.ts
  } else if (output_format == "tsibble") {
    as_tsibble(evolved.ts)
  }
  return(output)
}
