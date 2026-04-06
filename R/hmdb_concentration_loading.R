#' @title Open the HMDB concentration dataset 
#' @description
#'
#' `hmdb_concentration_loading` opens the HMDB concentration dataset as a tibble and return it.
#' 
#' @return a tibble
#' @export
hmdb_concentration_loading <- function() {
  hmdb_concentration_candidates <- massmatcher_inst_file_candidates(
    "concentration",
    "hmdb_concentrations_normal_condition.csv"
  )
  hmdb_concentration_path <- hmdb_concentration_candidates[file.exists(hmdb_concentration_candidates)][1]
  if (is.na(hmdb_concentration_path)) {
    stop(
      "Could not locate hmdb_concentrations_normal_condition.csv. Checked: ",
      paste(hmdb_concentration_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  cache_key <- paste0(
    "hmdb_concentration::",
    normalizePath(hmdb_concentration_path, winslash = "/", mustWork = FALSE)
  )
  cached <- massmatcher_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  message("HMDB concentration dataset is loading...")
  hmdb_concentration_normal <- readr::read_csv(
    hmdb_concentration_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  massmatcher_cache_set(cache_key, hmdb_concentration_normal)
}
