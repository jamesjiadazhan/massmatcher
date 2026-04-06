testthat::test_that("annotate_match_results enriches mass and isotopologue outputs", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("tibble")

  ref <- massmatcher:::metabolite_database_loading(database = "metorigindb") |>
    dplyr::filter(
      !is.na(.data$InChIKey),
      !is.na(.data$Subclass),
      .data$Subclass != "",
      .data$Charge_natural == 0,
      !is.na(.data$Mono_mass),
      !is.na(.data$Most_abundant_isotopologue_mass),
      !is.na(.data$Exact_mass),
      !is.na(.data$Exact_mass_most_abundant_isotopologue)
    ) |>
    dplyr::select(
      Mono_mass,
      Most_abundant_isotopologue_mass,
      Exact_mass,
      Exact_mass_most_abundant_isotopologue
    ) |>
    utils::head(1) |>
    dplyr::collect()

  if (nrow(ref) == 0) {
    testthat::skip("No neutral metorigindb reference molecule available.")
  }

  mono_mz <- massmatcher:::mass2mz_df_safe(ref$Exact_mass[[1]], adduct = "M+H")$mz[[1]]
  iso_mz <- massmatcher:::mass2mz_df_safe(ref$Exact_mass_most_abundant_isotopologue[[1]], adduct = "M+H")$mz[[1]]

  mono_out <- mass_match(
    unknown_feature = tibble::tibble(mz = mono_mz, time = 1),
    adduct = "M+H",
    database = "metorigindb"
  )
  iso_out <- isotopologue_mass_match(
    unknown_feature = tibble::tibble(mz = iso_mz, time = 1),
    adduct = "M+H",
    database = "metorigindb"
  )

  mono_out <- annotate_match_results(
    result_table = mono_out,
    database = "metorigindb",
    query_missing_classification = FALSE
  )
  iso_out <- annotate_match_results(
    result_table = iso_out,
    database = "metorigindb",
    query_missing_classification = FALSE
  )

  expected_cols <- c(
    "Classification_Kingdom",
    "Classification_Superclass",
    "Classification_Class",
    "Classification_Subclass",
    "Classification_Direct_parent",
    "Classification_Alternative_parent",
    "Concentration_average",
    "Concentration_units"
  )

  testthat::expect_true(nrow(mono_out) > 0)
  testthat::expect_true(nrow(iso_out) > 0)
  testthat::expect_true(all(expected_cols %in% colnames(mono_out)))
  testthat::expect_true(all(expected_cols %in% colnames(iso_out)))
  testthat::expect_true(any(!is.na(mono_out$Classification_Subclass)))
  testthat::expect_true(any(!is.na(iso_out$Classification_Subclass)))
})

testthat::test_that("annotate_mz_match_clustering_results enriches clustering outputs", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("preprocessCore")

  data("feature_table_exp_hilicpos", package = "massmatcher")

  feature_subset <- feature_table_exp_hilicpos[1:200, 1:10]
  feature_subset <- dplyr::as_tibble(feature_subset)
  colnames(feature_subset)[1:2] <- c("mz", "time")

  target_mz <- feature_subset$mz[1:3]
  synthetic_masses <- massmatcher:::mz2mass_df_safe(x = target_mz, adduct = "M+H")$adduct_mass

  synthetic_database <- tibble::tibble(
    Mono_mass = synthetic_masses,
    Most_abundant_isotopologue_mass = synthetic_masses,
    Exact_mass = synthetic_masses,
    Exact_mass_most_abundant_isotopologue = synthetic_masses,
    Name = paste0("Synthetic_", seq_along(synthetic_masses)),
    Formula = rep("C6H12O6", length(synthetic_masses)),
    HMDB_ID = rep(NA_character_, length(synthetic_masses)),
    KEGG_ID = rep(NA_character_, length(synthetic_masses)),
    InChIKey = sprintf("SYN%011d", seq_along(synthetic_masses)),
    Subclass = rep("Synthetic", length(synthetic_masses)),
    has_NH3 = rep(FALSE, length(synthetic_masses)),
    has_H2O = rep(FALSE, length(synthetic_masses)),
    n_NH3 = rep(0L, length(synthetic_masses)),
    n_H2O = rep(0L, length(synthetic_masses))
  )

  clustering_out <- mz_match_clustering(
    met_raw_wide = feature_subset,
    metabolite_database = synthetic_database,
    database = "metorigindb",
    mz_threshold = 10,
    All_Adduct = "M+H",
    imputation_method = "half_min",
    write_output = FALSE
  )

  out <- annotate_mz_match_clustering_results(
    clustering_output = clustering_out,
    database = "metorigindb",
    query_missing_classification = FALSE
  )

  testthat::expect_true(all(c("mz_only", "mz_only_isotope", "clustering_all", "clustering_main") %in% names(out)))

  for (nm in c("mz_only", "mz_only_isotope", "clustering_all", "clustering_main")) {
    tbl <- out[[nm]]
    testthat::expect_s3_class(tbl, "tbl_df")
    testthat::expect_true("Concentration_average" %in% colnames(tbl))
    testthat::expect_true("Classification_Subclass" %in% colnames(tbl))
  }
})
