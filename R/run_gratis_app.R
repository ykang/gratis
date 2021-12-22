#' Web Application to generate time series with controllable features.
#'
#' @return
#' NULL
#
# @examples
# # Not Run
# run_gratis_app <- function() {
#   appDir <- system.file("shiny", "gratis", package = "gratis")
#   if (appDir == "") {
#     stop("Could not find example directory. Try re-installing `gratis`.", call. = FALSE)
#   }
# 
#   shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)
# }

#' @rdname app
#' @export
app_gratis <- function() {
  appDir <- system.file("shiny", "app", package = "gratis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `gratis`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)
}
