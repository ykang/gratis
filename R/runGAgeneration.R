#' Title
#'
#' @return
#' @export
#'
#' @examples
runGAgeneration <- function() {
  appDir <- system.file("shiny-generation", "GAgeneration", package = "tsgeneration")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `tsgeneration`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)
}