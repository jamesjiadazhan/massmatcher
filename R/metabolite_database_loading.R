locate_database_parquet <- function(database = c("metorigindb", "pubchem")) {
  database <- match.arg(database)

  if (identical(database, "metorigindb")) {
    package_subdir <- "metorigindb"
    file_name <- "metorigindb_database_major_enriched.parquet"
    message_text <- "metorigindb dataset is loading..."
  } else {
    package_subdir <- "pubchem"
    file_name <- "pubchem_clean-0.parquet"
    message_text <- "PubChem dataset is loading..."
  }

  candidate_paths <- c(
    system.file(package_subdir, file_name, package = "massmatcher"),
    file.path(getwd(), "inst", package_subdir, file_name),
    file.path(getwd(), "massmatcher", "inst", package_subdir, file_name)
  )

  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])
  parquet_path <- candidate_paths[file.exists(candidate_paths)][1]

  if (is.na(parquet_path)) {
    stop(
      "Could not locate the ", database, " parquet file. Checked: ",
      paste(candidate_paths, collapse = ", "),
      call. = FALSE
    )
  }

  list(path = parquet_path, message = message_text)
}

metabolite_database_loading <- function(database = c("metorigindb", "pubchem")) {
  database <- match.arg(database)
  location <- locate_database_parquet(database)

  cache_key <- paste0(
    "metabolite_raw_dataset::",
    database,
    "::",
    normalizePath(location$path, winslash = "/", mustWork = FALSE)
  )
  dataset <- massmatcher_cache_get(cache_key)
  if (is.null(dataset)) {
    message(location$message)
    dataset <- arrow::open_dataset(location$path)
    massmatcher_cache_set(cache_key, dataset)
  }

  if (identical(database, "metorigindb")) {
    dataset |>
      dplyr::transmute(
        Name = .data$Compound_name,
        Formula = .data$Formula,
        Mono_mass = .data$Monoisotopic_weight,
        Most_abundant_isotopologue_mass = .data$Most_abundant_isotopologue_mass,
        InChIKey = .data$INCHIKEY_ID,
        SMILES = .data$SMILES_ID,
        Charge_natural = .data$Charge_natural,
        unique_identifier = .data$unique_identifier,
        HMDB_ID = .data$HMDB_ID,
        KEGG_ID = .data$KEGG_ID,
        CAS_ID = .data$CAS_ID,
        INCHI_ID = .data$INCHI_ID,
        KEGG_DRUG_ID = .data$KEGG_DRUG_ID,
        CHEMSPIDER_ID = .data$CHEMSPIDER_ID,
        DRUGBANK_ID = .data$DRUGBANK_ID,
        FOODB_ID = .data$FOODB_ID,
        PUBCHEM_COMPOUND_ID = .data$PUBCHEM_COMPOUND_ID,
        PUBCHEM_SUBSTANCE_ID = .data$PUBCHEM_SUBSTANCE_ID,
        CHEBI_ID = .data$CHEBI_ID,
        CHEMBL_ID = .data$CHEMBL_ID,
        PDB_CCD_ID = .data$PDB_CCD_ID,
        `3DMET_ID` = .data$`3DMET_ID`,
        NIKKAJI_ID = .data$NIKKAJI_ID,
        KNAPSACK_ID = .data$KNAPSACK_ID,
        LIPIDMAPS_ID = .data$LIPIDMAPS_ID,
        LIPIDBANK_ID = .data$LIPIDBANK_ID,
        BIOCYC_ID = .data$BIOCYC_ID,
        BIGG_ID = .data$BIGG_ID,
        BIGG_IDENTIFIER_ID = .data$BIGG_IDENTIFIER_ID,
        WIKIPEDIA_ID = .data$WIKIPEDIA_ID,
        METLIN_ID = .data$METLIN_ID,
        T3DB_ID = .data$T3DB_ID,
        REACTOME_ID = .data$REACTOME_ID,
        MODELSEED_ID = .data$MODELSEED_ID,
        MIMEDB_ID = .data$MIMEDB_ID,
        LOTUS_ID = .data$LOTUS_ID,
        Subclass = .data$Subclass,
        has_NH3 = .data$has_NH3,
        n_NH3 = .data$n_NH3,
        has_H2O = .data$has_H2O,
        n_H2O = .data$n_H2O,
        Mass_diff = .data$Mass_diff
      )
  } else {
    dataset |>
      dplyr::transmute(
        Name = .data$Name,
        Formula = .data$Formula,
        Mono_mass = .data$MonoMass,
        Most_abundant_isotopologue_mass = .data$most_abundant_isotopologue_mass,
        InChIKey = .data$InChIKey,
        InChI = .data$InChI,
        SMILES = .data$SMILES,
        Charge_natural = .data$Charge_natural,
        CID = .data$CID,
        HMDB_ID = .data$HMDB_ID,
        KEGG_ID = .data$KEGG_ID,
        LIPID_MAPS_ID = .data$LIPID_MAPS_ID,
        ChEBI_ID = .data$ChEBI_ID,
        BioCyc_ID = .data$BioCyc_ID,
        DrugBank_ID = .data$DrugBank_ID,
        Subclass = NA_character_,
        has_NH3 = FALSE,
        n_NH3 = 0L,
        has_H2O = FALSE,
        n_H2O = 0L,
        Mass_diff = NA_real_
      )
  }
}

adduct_definition_loading <- function() {
  adduct_env <- new.env(parent = emptyenv())
  utils::data("adduct_definition", package = "MetaboCoreUtilsAdduct", envir = adduct_env)

  if (!exists("adduct_definition", envir = adduct_env, inherits = FALSE)) {
    stop("Could not load adduct_definition from MetaboCoreUtilsAdduct.", call. = FALSE)
  }

  get("adduct_definition", envir = adduct_env, inherits = FALSE)
}

mass2mz_df_safe <- function(mass, adduct = "M+H") {
  adduct_definition <- adduct_definition_loading()
  mass_adduct <- data.frame(mass = mass, adduct = adduct, stringsAsFactors = FALSE)

  mass_adduct_adduct_definition <- dplyr::left_join(
    mass_adduct,
    adduct_definition,
    by = c("adduct" = "name")
  )

  mass_adduct_adduct_definition$mz <- (
    mass_adduct_adduct_definition$mass * mass_adduct_adduct_definition$mass_multi
  ) + mass_adduct_adduct_definition$mass_add

  mass_adduct_adduct_definition[, c("mass", "adduct", "mz")]
}

mz2mass_df_safe <- function(x, adduct = "M+H") {
  adduct_definition <- adduct_definition_loading() |>
    dplyr::filter(.data$name %in% adduct)

  comb <- expand.grid(
    x = x,
    name = adduct_definition$name,
    stringsAsFactors = FALSE
  )
  comb <- merge(comb, adduct_definition, by = "name")
  comb$adduct_mass <- (comb$x - comb$mass_add) / comb$mass_multi
  comb <- comb[, c("x", "name", "adduct_mass")]
  names(comb) <- c("mz", "adduct", "adduct_mass")
  comb
}
