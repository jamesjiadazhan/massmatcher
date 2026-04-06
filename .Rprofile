profile_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NULL
)
profile_dir <- if (is.null(profile_path)) getwd() else dirname(profile_path)

renv_activate <- file.path(profile_dir, "renv", "activate.R")
if (file.exists(renv_activate)) {
  source(renv_activate)
}
