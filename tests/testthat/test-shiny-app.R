testthat::test_that("bundled shiny app file is present", {
  app_file <- system.file("shiny", "massmatcher-app", "app.R", package = "massmatcher")
  testthat::expect_true(nzchar(app_file))
  testthat::expect_true(file.exists(app_file))
})

testthat::test_that("clustering demo controls are present in shiny app", {
  app_file <- system.file("shiny", "massmatcher-app", "app.R", package = "massmatcher")
  app_code <- readLines(app_file, warn = FALSE)

  testthat::expect_true(any(grepl("cluster_feature_source_ui", app_code, fixed = TRUE)))
  testthat::expect_true(any(grepl("download_cluster_demo", app_code, fixed = TRUE)))
  testthat::expect_true(any(grepl("default_cluster_adducts", app_code, fixed = TRUE)))
  testthat::expect_true(any(grepl("inputId = \"cluster_adducts\"", app_code, fixed = TRUE)))
  testthat::expect_false(any(grepl("cluster_biospecimen", app_code, fixed = TRUE)))
  testthat::expect_true(any(grepl("inputId = \"meta_biospecimen\"", app_code, fixed = TRUE)))
  testthat::expect_true(any(grepl("metadata_download_ui", app_code, fixed = TRUE)))
  testthat::expect_false(any(grepl("meta_query_missing", app_code, fixed = TRUE)))
  testthat::expect_false(any(grepl("meta_match_first_block", app_code, fixed = TRUE)))
})

testthat::test_that("bundled shiny logo file is present", {
  logo_file <- system.file("shiny", "massmatcher-app", "www", "massmatcher-logo.svg", package = "massmatcher")
  testthat::expect_true(nzchar(logo_file))
  testthat::expect_true(file.exists(logo_file))
})
