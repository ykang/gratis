#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
nsdiffs1 <- function(x){
  c(nsdiffs=ifelse(frequency(x)==1L, -1, forecast::nsdiffs(x)))
}