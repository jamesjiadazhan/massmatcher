testthat::test_that("mz_match_clustering works with feature_table_exp_hilicpos example data", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("preprocessCore")

  data("feature_table_exp_hilicpos", package = "massmatcher")

  feature_subset <- feature_table_exp_hilicpos[1:250, 1:10]
  feature_subset <- dplyr::as_tibble(feature_subset)
  colnames(feature_subset)[1:2] <- c("mz", "time")

  target_mz <- feature_subset$mz[1:3]
  synthetic_masses <- massmatcher:::mz2mass_df_safe(x = target_mz, adduct = "M+H")$adduct_mass

  synthetic_database <- tibble::tibble(
    Mono_mass = synthetic_masses,
    Most_abundant_isotopologue_mass = synthetic_masses,
    Name = paste0("Synthetic_", seq_along(synthetic_masses)),
    Formula = rep("C6H12O6", length(synthetic_masses)),
    HMDB_ID = rep(NA_character_, length(synthetic_masses)),
    KEGG_ID = rep(NA_character_, length(synthetic_masses)),
    InChIKey = sprintf("SYNTHETICKEY%02d-ABCDEFGHIJ-A", seq_along(synthetic_masses)),
    Subclass = rep("Synthetic", length(synthetic_masses)),
    has_NH3 = rep(FALSE, length(synthetic_masses)),
    has_H2O = rep(FALSE, length(synthetic_masses)),
    n_NH3 = rep(0L, length(synthetic_masses)),
    n_H2O = rep(0L, length(synthetic_masses))
  )

  result <- mz_match_clustering(
    met_raw_wide = feature_subset,
    metabolite_database = synthetic_database,
    database = "metorigindb",
    mz_threshold = 10,
    All_Adduct = "M+H",
    imputation_method = "half_min",
    write_output = FALSE
  )

  testthat::expect_true(is.list(result))
  testthat::expect_true(all(c("mz_only", "mz_only_isotope", "clustering_all", "clustering_main") %in% names(result)))
  testthat::expect_s3_class(result$mz_only, "tbl_df")
  testthat::expect_true(nrow(result$mz_only) > 0)
  testthat::expect_true(all(c("unique_identifier", "InChIKey", "SMILES", "KEGG_ID", "HMDB_ID") %in% colnames(result$mz_only)))
  testthat::expect_true(all(c("unique_identifier", "InChIKey", "SMILES", "KEGG_ID", "HMDB_ID") %in% colnames(result$mz_only_isotope)))
  testthat::expect_true(all(c("unique_identifier", "InChIKey", "SMILES", "KEGG_ID", "HMDB_ID") %in% colnames(result$clustering_all)))
  testthat::expect_true(all(c("unique_identifier", "InChIKey", "SMILES", "KEGG_ID", "HMDB_ID") %in% colnames(result$clustering_main)))
  testthat::expect_false("CID" %in% colnames(result$mz_only))
  testthat::expect_false("CID" %in% colnames(result$mz_only_isotope))
  testthat::expect_false("CID" %in% colnames(result$clustering_all))
  testthat::expect_false("CID" %in% colnames(result$clustering_main))
  testthat::expect_lt(
    match("theoretical_mz", colnames(result$mz_only)),
    match("mz_matching_ppm", colnames(result$mz_only))
  )
  testthat::expect_lt(
    match("theoretical_mz", colnames(result$mz_only_isotope)),
    match("mz_matching_ppm", colnames(result$mz_only_isotope))
  )

  intensity_columns <- colnames(feature_subset)[-c(1, 2)]
  testthat::expect_false(any(intensity_columns %in% colnames(result$mz_only)))
  testthat::expect_false(any(intensity_columns %in% colnames(result$mz_only_isotope)))
  testthat::expect_false(any(intensity_columns %in% colnames(result$clustering_all)))
  testthat::expect_false(any(intensity_columns %in% colnames(result$clustering_main)))
  testthat::expect_true(all(result$mz_only$InChIKey %in% synthetic_database$InChIKey))
  testthat::expect_true(any(nchar(result$mz_only$InChIKey) > 14))

  if (all(c("mz_time_isotope_annotated", "mz_time_annotated") %in% colnames(result$clustering_all))) {
    testthat::expect_false(any(
      result$clustering_all$mz_time_isotope_annotated == result$clustering_all$mz_time_annotated,
      na.rm = TRUE
    ))
  }

  if (all(c("mz_time_isotope_annotated", "mz_time_annotated") %in% colnames(result$clustering_main))) {
    testthat::expect_false(any(
      result$clustering_main$mz_time_isotope_annotated == result$clustering_main$mz_time_annotated,
      na.rm = TRUE
    ))
  }
})
