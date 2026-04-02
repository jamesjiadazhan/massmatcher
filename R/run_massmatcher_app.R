#' @title Launch the massmatcher Shiny app
#' @description
#' Starts the bundled Shiny app that provides three major modules:
#' Search, Clustering, and Additional Meta Data enrichment.
#'
#' @param host Host to bind the Shiny app. Default is `"127.0.0.1"`.
#' @param port Port for the Shiny app. Default is `NULL` (auto-select).
#' @param launch.browser Logical; whether to open a browser automatically.
#'   Default is `TRUE`.
#'
#' @return The result of `shiny::runApp()`.
#' @export
run_massmatcher_app <- function(host = "127.0.0.1", port = NULL, launch.browser = TRUE) {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Please install the `shiny` package to run the app.", call. = FALSE)
    }
    if (!requireNamespace("bslib", quietly = TRUE)) {
        stop("Please install the `bslib` package to run the app.", call. = FALSE)
    }
    if (!requireNamespace("DT", quietly = TRUE)) {
        stop("Please install the `DT` package to run the app.", call. = FALSE)
    }

    app_dir <- system.file("shiny", "massmatcher-app", package = "massmatcher")
    if (identical(app_dir, "")) {
        stop("Shiny app directory not found in the installed package.", call. = FALSE)
    }

    shiny::runApp(
        appDir = app_dir,
        host = host,
        port = port,
        launch.browser = launch.browser
    )
}
