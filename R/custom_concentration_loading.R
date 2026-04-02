#' @title Open the custom concentration dataset 
#' @description
#'
#' `custom_concentration_loading` opens the custom concentration dataset (provided by the Jones/Go Lab, Department of Medicine, Emory University) as a tibble and return it. User can add their own concentration dataset in the same format and load it with this function.
#' 
#' @return a tibble
#' @export
custom_concentration_loading <- function() {
  custom_concentration_path <- system.file("concentration", "custom_concentrations.csv", package = "massmatcher")
  if (custom_concentration_path == "") {
    stop(
      "Could not locate custom_concentrations.csv in the installed package.",
      call. = FALSE
    )
  }

  cache_key <- paste0(
    "custom_concentration::",
    normalizePath(custom_concentration_path, winslash = "/", mustWork = FALSE)
  )
  cached <- massmatcher_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  message("Custom concentration dataset is loading...")
  custom_concentration <- readr::read_csv(
    custom_concentration_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  massmatcher_cache_set(cache_key, custom_concentration)
}
