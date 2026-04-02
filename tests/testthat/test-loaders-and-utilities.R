testthat::test_that("concentration loaders return non-empty tables", {
  testthat::skip_if_not_installed("readr")
  testthat::skip_if_not_installed("tibble")

  hmdb <- hmdb_concentration_loading()
  custom <- custom_concentration_loading()

  testthat::expect_s3_class(hmdb, "tbl_df")
  testthat::expect_true(nrow(hmdb) > 0)
  testthat::expect_true(all(c("HMDBID", "Concentration_average", "Concentration_units") %in% colnames(hmdb)))

  testthat::expect_s3_class(custom, "tbl_df")
  testthat::expect_true(nrow(custom) > 0)
})

testthat::test_that("pubchem_classyfire_loading returns expected dataset schema", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")

  classyfire_cols <- pubchem_classyfire_loading() |>
    utils::head(0) |>
    dplyr::collect() |>
    colnames()

  testthat::expect_true("InChIKey" %in% classyfire_cols)
  testthat::expect_true("subclass" %in% classyfire_cols)
})

testthat::test_that("find.Overlapping.mzs finds overlaps with and without time threshold", {
  mz_overlap <- find.Overlapping.mzs(
    dataA = data.frame(mz = c(100, 200)),
    dataB = data.frame(mz = c(100.0002, 300)),
    mz.thresh = 5
  )
  testthat::expect_true(nrow(mz_overlap) == 1)
  testthat::expect_equal(as.integer(mz_overlap$index.A[[1]]), 1L)
  testthat::expect_equal(as.integer(mz_overlap$index.B[[1]]), 1L)

  no_time_overlap <- find.Overlapping.mzs(
    dataA = data.frame(mz = 100, time = 10),
    dataB = data.frame(mz = 100.0002, time = 50),
    mz.thresh = 5,
    time.thresh = 30
  )
  testthat::expect_true(nrow(no_time_overlap) == 0)

  with_time_overlap <- find.Overlapping.mzs(
    dataA = data.frame(mz = 100, time = 10),
    dataB = data.frame(mz = 100.0002, time = 20),
    mz.thresh = 5,
    time.thresh = 30
  )
  testthat::expect_true(nrow(with_time_overlap) == 1)
})

testthat::test_that("find.Overlapping.mzs can recover numeric m/z columns after coercion", {
  data_a <- data.frame(
    mz = c("100.0000", "200.0000"),
    stringsAsFactors = FALSE
  )
  data_b <- data.frame(
    label = c("candidate_1", "candidate_2"),
    mz = c("100.0003", "not-a-number"),
    stringsAsFactors = FALSE
  )

  overlap <- find.Overlapping.mzs(
    dataA = data_a,
    dataB = data_b,
    mz.thresh = 5
  )

  testthat::expect_equal(nrow(overlap), 1)
  testthat::expect_equal(as.integer(overlap$index.A[[1]]), 1L)
  testthat::expect_equal(as.integer(overlap$index.B[[1]]), 1L)
})

testthat::test_that("target_mass_search recovers expected feature match", {
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("tibble")

  target_mass <- 180.063388
  target_mz <- massmatcher:::mass2mz_df_safe(mass = target_mass, adduct = "M+H")$mz[[1]]
  feature_table <- tibble::tibble(
    mz = c(target_mz, target_mz + 1),
    time = c(100, 120),
    intensity = c(10000, 5000)
  )

  out <- target_mass_search(
    mass = target_mass,
    feature_table = feature_table,
    adduct = "M+H",
    mz_ppm = 10
  )

  testthat::expect_s3_class(out, "tbl_df")
  testthat::expect_true(nrow(out) >= 1)
  testthat::expect_true(min(out$mz_ppm, na.rm = TRUE) <= 0.01)
})
