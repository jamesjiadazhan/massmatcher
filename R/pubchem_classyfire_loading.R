#' Open the preprocessed PubChem ClassyFire dataset 
#'
#' `pubchem_classyfire_loading` opens the PubChem ClassyFire dataset (with ClassyFire classification) preprocessed as an Arrow Dataset, so you can use dplyr verbs without loading the whole file to avoid memory burden.
#'
#' @return an `arrow::Dataset`
#' @export 
pubchem_classyfire_loading <- function() {
  parquet_candidates <- massmatcher_inst_file_candidates("classyfire", "pubchem_classyfire-0.parquet")
  parquet_path <- parquet_candidates[file.exists(parquet_candidates)][1]
  if (is.na(parquet_path)) {
    stop(
      "Could not locate pubchem_classyfire-0.parquet. Checked: ",
      paste(parquet_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  cache_key <- paste0(
    "pubchem_classyfire_dataset::",
    normalizePath(parquet_path, winslash = "/", mustWork = FALSE)
  )
  cached <- massmatcher_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  message("Pubchem ClassyFire dataset is loading...")
  dataset <- arrow::open_dataset(parquet_path)
  massmatcher_cache_set(cache_key, dataset)
}
