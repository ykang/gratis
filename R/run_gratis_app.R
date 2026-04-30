#' Launch the gratis Shiny application
#'
#' Launch a local Shiny application for generating time series with
#' controllable features.
#'
#' @return Called for its side effect of launching a Shiny application. Returns
#'   the value from \code{\link[shiny]{runApp}()}.
#
#' @examples
#' \dontrun{
#' app_gratis()
#' }
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
