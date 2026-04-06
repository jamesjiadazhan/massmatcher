massmatcher_bundle_root <- function() {
  root <- getOption("massmatcher.bundle_root", default = Sys.getenv("MASSMATCHER_BUNDLE_ROOT", unset = ""))
  if (!is.character(root) || length(root) == 0L) {
    return(NULL)
  }

  root <- trimws(root[[1]])
  if (!nzchar(root) || !dir.exists(root)) {
    return(NULL)
  }

  normalizePath(root, winslash = "/", mustWork = FALSE)
}

massmatcher_inst_file_candidates <- function(subdir, file_name) {
  bundle_root <- massmatcher_bundle_root()

  candidates <- c(
    system.file(subdir, file_name, package = "massmatcher"),
    if (!is.null(bundle_root)) file.path(bundle_root, "inst", subdir, file_name),
    file.path(getwd(), "inst", subdir, file_name),
    file.path(getwd(), "massmatcher", "inst", subdir, file_name)
  )

  unique(candidates[nzchar(candidates)])
}

massmatcher_data_file_candidates <- function(file_name) {
  bundle_root <- massmatcher_bundle_root()

  candidates <- c(
    system.file("data", file_name, package = "massmatcher"),
    if (!is.null(bundle_root)) file.path(bundle_root, "data", file_name),
    file.path(getwd(), "data", file_name),
    file.path(getwd(), "massmatcher", "data", file_name)
  )

  unique(candidates[nzchar(candidates)])
}
