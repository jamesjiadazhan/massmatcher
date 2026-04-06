testthat::test_that("database loaders support metorigindb and pubchem", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")

  metorigindb_cols <- pubchem_loading() |>
    utils::head(0) |>
    dplyr::collect() |>
    colnames()
  pubchem_cols <- pubchem_loading(database = "pubchem") |>
    utils::head(0) |>
    dplyr::collect() |>
    colnames()

  testthat::expect_true("Compound_name" %in% metorigindb_cols)
  testthat::expect_true("Name" %in% pubchem_cols)
  testthat::expect_true(all(c("Exact_mass", "Exact_mass_most_abundant_isotopologue") %in% colnames(
    massmatcher:::metabolite_database_loading(database = "metorigindb") |>
      utils::head(0) |>
      dplyr::collect()
  )))
  testthat::expect_true(all(c("Exact_mass", "Exact_mass_most_abundant_isotopologue") %in% colnames(
    massmatcher:::metabolite_database_loading(database = "pubchem") |>
      utils::head(0) |>
      dplyr::collect()
  )))
})

testthat::test_that("database-backed functions default to metorigindb", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("tibble")

  reference <- massmatcher:::metabolite_database_loading(database = "metorigindb") |>
    dplyr::filter(
      !is.na(.data$InChIKey),
      .data$Charge_natural == 0,
      !is.na(.data$Mono_mass),
      !is.na(.data$Exact_mass)
    ) |>
    dplyr::select(InChIKey, Mono_mass, Exact_mass) |>
    utils::head(1) |>
    dplyr::collect()

  if (nrow(reference) == 0) {
    testthat::skip("No neutral metorigindb reference molecule available.")
  }

  mono_mz <- massmatcher:::mass2mz_df_safe(
    mass = reference$Exact_mass[[1]],
    adduct = "M+H"
  )$mz[[1]]
  ref_key <- toupper(trimws(reference$InChIKey[[1]]))

  default_class <- get_chemical_classification(
    unique_inchikey = reference$InChIKey,
    query_missing = FALSE
  )
  explicit_class <- get_chemical_classification(
    unique_inchikey = reference$InChIKey,
    query_missing = FALSE,
    database = "metorigindb"
  )

  testthat::expect_equal(default_class$InChIKey, explicit_class$InChIKey)
  testthat::expect_equal(default_class$Subclass, explicit_class$Subclass)
})

testthat::test_that("mass and isotopologue matching functions work for both databases", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("tibble")

  for (database_name in c("metorigindb", "pubchem")) {
    reference <- massmatcher:::metabolite_database_loading(database = database_name) |>
      dplyr::filter(
        !is.na(.data$InChIKey),
        .data$Charge_natural == 0,
        !is.na(.data$Mono_mass),
        !is.na(.data$Most_abundant_isotopologue_mass),
        !is.na(.data$Exact_mass),
        !is.na(.data$Exact_mass_most_abundant_isotopologue)
      ) |>
      dplyr::select(
        InChIKey,
        Mono_mass,
        Most_abundant_isotopologue_mass,
        Exact_mass,
        Exact_mass_most_abundant_isotopologue
      ) |>
      utils::head(1) |>
      dplyr::collect()

    if (nrow(reference) == 0) {
      testthat::skip(paste0("No neutral reference molecule available in ", database_name, "."))
    }

    ref_key <- toupper(trimws(reference$InChIKey[[1]]))
    mono_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Exact_mass[[1]],
      adduct = "M+H"
    )$mz[[1]]
    isotopologue_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Exact_mass_most_abundant_isotopologue[[1]],
      adduct = "M+H"
    )$mz[[1]]

    mass_matches <- mass_match(
      unknown_feature = tibble::tibble(mz = mono_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(mass_matches) > 0)
    testthat::expect_true(any(toupper(trimws(mass_matches$InChIKey)) == ref_key))
    testthat::expect_false(any(c("Exact_mass", "Exact_mass_most_abundant_isotopologue") %in% colnames(mass_matches)))

    isotopologue_matches <- isotopologue_mass_match(
      unknown_feature = tibble::tibble(mz = isotopologue_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(isotopologue_matches) > 0)
    testthat::expect_true(any(toupper(trimws(isotopologue_matches$InChIKey)) == ref_key))
    testthat::expect_false(any(c("Exact_mass", "Exact_mass_most_abundant_isotopologue") %in% colnames(isotopologue_matches)))
  }
})

testthat::test_that("package-level matching defaults align with explicit broad prefiltering", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MetaboCoreUtilsAdduct")
  testthat::skip_if_not_installed("tibble")

  for (database_name in c("metorigindb", "pubchem")) {
    reference <- massmatcher:::metabolite_database_loading(database = database_name) |>
      dplyr::filter(
        !is.na(.data$InChIKey),
        .data$Charge_natural == 0,
        !is.na(.data$Exact_mass),
        !is.na(.data$Exact_mass_most_abundant_isotopologue)
      ) |>
      dplyr::select(InChIKey, Exact_mass, Exact_mass_most_abundant_isotopologue) |>
      utils::head(1) |>
      dplyr::collect()

    if (nrow(reference) == 0) {
      testthat::skip(paste0("No neutral reference molecule available in ", database_name, "."))
    }

    ref_key <- toupper(trimws(reference$InChIKey[[1]]))
    mono_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Exact_mass[[1]],
      adduct = "M+H"
    )$mz[[1]]
    isotopologue_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Exact_mass_most_abundant_isotopologue[[1]],
      adduct = "M+H"
    )$mz[[1]]

    mono_prefiltered <- prefilter_metabolite_database_by_mz_range(
      mz_values = mono_mz,
      adducts = "M+H",
      ppm_threshold = 10,
      database = database_name
    )
    mono_default <- mass_match(
      unknown_feature = tibble::tibble(mz = mono_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    mono_explicit <- mass_match(
      unknown_feature = tibble::tibble(mz = mono_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name,
      metabolite_database = mono_prefiltered
    )

    iso_prefiltered <- prefilter_metabolite_database_by_mz_range(
      mz_values = isotopologue_mz,
      adducts = "M+H",
      ppm_threshold = 10,
      database = database_name
    )
    iso_default <- isotopologue_mass_match(
      unknown_feature = tibble::tibble(mz = isotopologue_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    iso_explicit <- isotopologue_mass_match(
      unknown_feature = tibble::tibble(mz = isotopologue_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name,
      metabolite_database = iso_prefiltered
    )

    testthat::expect_true(any(toupper(trimws(mono_default$InChIKey)) == ref_key))
    testthat::expect_true(any(toupper(trimws(mono_explicit$InChIKey)) == ref_key))
    testthat::expect_true(any(toupper(trimws(iso_default$InChIKey)) == ref_key))
    testthat::expect_true(any(toupper(trimws(iso_explicit$InChIKey)) == ref_key))
  }
})

testthat::test_that("classification cache supports metorigindb and pubchem", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("classyfireR")

  metorigindb_ref <- massmatcher:::metabolite_database_loading(database = "metorigindb") |>
    dplyr::filter(!is.na(.data$InChIKey), !is.na(.data$Subclass), .data$Subclass != "") |>
    dplyr::select(InChIKey, Subclass) |>
    utils::head(1) |>
    dplyr::collect()

  if (nrow(metorigindb_ref) == 0) {
    testthat::skip("No metorigindb row with non-missing Subclass available.")
  }

  metorigindb_class <- get_chemical_classification(
    unique_inchikey = metorigindb_ref$InChIKey,
    query_missing = FALSE,
    database = "metorigindb"
  )
  metorigindb_class_pubchem_arg <- get_chemical_classification(
    unique_inchikey = metorigindb_ref$InChIKey,
    query_missing = FALSE,
    database = "pubchem"
  )
  testthat::expect_true(nrow(metorigindb_class) == 1)
  testthat::expect_equal(metorigindb_class$InChIKey[[1]], metorigindb_ref$InChIKey[[1]])
  testthat::expect_equal(metorigindb_class, metorigindb_class_pubchem_arg)

  pubchem_ref <- pubchem_classyfire_loading() |>
    dplyr::filter(!is.na(.data$InChIKey)) |>
    dplyr::select(InChIKey) |>
    utils::head(1) |>
    dplyr::collect()

  if (nrow(pubchem_ref) == 0) {
    testthat::skip("No pubchem classyfire row available.")
  }

  pubchem_class <- get_chemical_classification(
    unique_inchikey = pubchem_ref$InChIKey,
    query_missing = FALSE,
    database = "pubchem"
  )
  testthat::expect_true(pubchem_ref$InChIKey[[1]] %in% pubchem_class$InChIKey)
  testthat::expect_true(all(c("Kingdom", "Superclass", "Class", "Subclass") %in% colnames(pubchem_class)))

  pubchem_fake_key <- "ZZZZZZZZZZZZZZ-FAKEFAKEFA-N"
  pubchem_fake_class <- get_chemical_classification(
    unique_inchikey = pubchem_fake_key,
    query_missing = FALSE,
    database = "pubchem"
  )
  testthat::expect_true(pubchem_fake_key %in% pubchem_fake_class$InChIKey)
  testthat::expect_true(all(is.na(pubchem_fake_class$Subclass)))
})
