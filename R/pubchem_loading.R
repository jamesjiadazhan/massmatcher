#' Open a preprocessed metabolite dataset
#'
#' `pubchem_loading` opens the selected metabolite dataset as an Arrow Dataset,
#' so you can use dplyr verbs without loading the whole file into memory.
#'
#' @param database Database used for metabolite lookup. Choose either
#'   `"metorigindb"` or `"pubchem"`. Default is `"metorigindb"`.
#' @return an `arrow::Dataset`
#' @export
pubchem_loading <- function(database = c("metorigindb", "pubchem")) {
  database <- match.arg(database)
  location <- locate_database_parquet(database)

  cache_key <- paste0(
    "metabolite_raw_dataset::",
    database,
    "::",
    normalizePath(location$path, winslash = "/", mustWork = FALSE)
  )
  cached <- massmatcher_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  message(location$message)
  dataset <- arrow::open_dataset(location$path)
  massmatcher_cache_set(cache_key, dataset)
}
