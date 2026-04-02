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
      !is.na(.data$Mono_mass)
    ) |>
    dplyr::select(InChIKey, Mono_mass) |>
    utils::head(1) |>
    dplyr::collect()

  if (nrow(reference) == 0) {
    testthat::skip("No neutral metorigindb reference molecule available.")
  }

  mono_mz <- massmatcher:::mass2mz_df_safe(
    mass = reference$Mono_mass[[1]],
    adduct = "M+H"
  )$mz[[1]]
  ref_key <- toupper(trimws(reference$InChIKey[[1]]))

  default_candidates <- mass_filter(mz = mono_mz, mz_ppm = 10, adduct = "M+H")
  explicit_candidates <- mass_filter(
    mz = mono_mz,
    mz_ppm = 10,
    adduct = "M+H",
    database = "metorigindb"
  )

  testthat::expect_true(any(toupper(trimws(default_candidates$InChIKey)) == ref_key))
  testthat::expect_true(any(toupper(trimws(explicit_candidates$InChIKey)) == ref_key))

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
        !is.na(.data$Most_abundant_isotopologue_mass)
      ) |>
      dplyr::select(InChIKey, Mono_mass, Most_abundant_isotopologue_mass) |>
      utils::head(1) |>
      dplyr::collect()

    if (nrow(reference) == 0) {
      testthat::skip(paste0("No neutral reference molecule available in ", database_name, "."))
    }

    ref_key <- toupper(trimws(reference$InChIKey[[1]]))
    mono_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Mono_mass[[1]],
      adduct = "M+H"
    )$mz[[1]]
    isotopologue_mz <- massmatcher:::mass2mz_df_safe(
      mass = reference$Most_abundant_isotopologue_mass[[1]],
      adduct = "M+H"
    )$mz[[1]]

    mass_candidates <- mass_filter(
      mz = mono_mz,
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(mass_candidates) > 0)
    testthat::expect_true(any(toupper(trimws(mass_candidates$InChIKey)) == ref_key))

    mass_matches <- mass_match(
      unknown_feature = tibble::tibble(mz = mono_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(mass_matches) > 0)
    testthat::expect_true(any(toupper(trimws(mass_matches$InChIKey)) == ref_key))

    isotopologue_candidates <- isotopologue_mass_filter(
      mz = isotopologue_mz,
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(isotopologue_candidates) > 0)
    testthat::expect_true(any(toupper(trimws(isotopologue_candidates$InChIKey)) == ref_key))

    isotopologue_matches <- isotopologue_mass_match(
      unknown_feature = tibble::tibble(mz = isotopologue_mz, time = 1),
      mz_ppm = 10,
      adduct = "M+H",
      database = database_name
    )
    testthat::expect_true(nrow(isotopologue_matches) > 0)
    testthat::expect_true(any(toupper(trimws(isotopologue_matches$InChIKey)) == ref_key))
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
