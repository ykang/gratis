#' Title
#'
#' @param n 
#' @param selected.features 
#' @param target 
#' @param seasonal 
#' @param ts.length 
#' @param freq 
#'
#' @return
#' @export
#'
#' @examples
generate_ts <- function(n, ts.length, freq, seasonal, features, selected.features, target){
  ga_min <-
    if (seasonal == 0){
      c(rep(0, 10))
    }else if (seasonal == 1){
      c(rep(0, 17))
    }else{
      c(rep(0, 35))
    }
  ga_max <- 
    if (seasonal == 0){
      c(rep(1, 10))
    }else if (seasonal == 1){
      c(rep(1, 17))
    }else{
      c(rep(1, 35))
    }
  evolved.ts <- c()
  while(ifelse(is.null(dim(evolved.ts)), 0 < 1, dim(evolved.ts)[2] < n)){
    GA <- myGA(type = "real-valued", fitness = tsgeneration::fitness, features = features, seasonal = seasonal, 
               ts.length, freq, target, 3, selected.features,
               n = ts.length,
               min = ga_min, 
               max = ga_max, 
               parallel = TRUE, popSize = 30, maxiter = 100, 
               pmutation = 0.3, pcrossover = 0.8, maxFitness = -0.05, 
               run = 30, keepBest = TRUE, monitor = GA::gaMonitor2)
    evolved.ts.new <- 
      unique(do.call(cbind,
                     eval(parse(text = paste( "list(", paste("GA@bestSol[[GA@iter - ", 0:(GA@run-1), ']]', sep = '', collapse = ','), ')')))), MARGIN = 2)
    evolved.ts <- cbind(evolved.ts, evolved.ts.new)
  }
  return(evolved.ts)
}