#' Web Application to generate time series with controllable features.
#'
#' @return
#' NULL
#' @export
#'
#' @examples
#' # Not Run
run_tsgeneration_app <- function() {
  appDir <- system.file("shiny", "tsgeneration", package = "tsgeneration")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `tsgeneration`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)
}
