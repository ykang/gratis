#' Title
#'
#' @param x
#'
#' @return
#' NA
#' @export
#'
#' @examples
#' # Not Run
nsdiffs1 <- function(x){
  c(nsdiffs=ifelse(frequency(x)==1L, -1, forecast::nsdiffs(x)))
}
