.massmatcher_dataset_cache <- new.env(parent = emptyenv())

massmatcher_cache_get <- function(key) {
  if (!exists(key, envir = .massmatcher_dataset_cache, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = .massmatcher_dataset_cache, inherits = FALSE)
}

massmatcher_cache_set <- function(key, value) {
  assign(key, value, envir = .massmatcher_dataset_cache)
  value
}
