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
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the gratis app.", call. = FALSE)
  }
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required to run the gratis app.", call. = FALSE)
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("Package 'rlang' is required to run the gratis app.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to run the gratis app.", call. = FALSE)
  }
  if (!requireNamespace("tsfeatures", quietly = TRUE)) {
    stop("Package 'tsfeatures' is required to run the gratis app.", call. = FALSE)
  }

  appDir <- system.file("shiny", "app", package = "gratis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `gratis`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)
}
