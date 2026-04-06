#' @title Run m/z Matching and Adduct/Isotope Clustering
#' @description
#' `mz_match_clustering()` reuses the preprocessing, m/z matching, and
#' adduct/isotope clustering workflow inside the `massmatcher` package.
#' It intentionally omits class-file handling, retention-time mapping,
#' and precursor-product/transporter logic.
#'
#' The function writes four CSV files by default:
#' `positive_metabolite_annotated_10ppm_mzonly.csv`,
#' `positive_metabolite_annotated_10ppm_mzonly_isotope.csv`,
#' `clustering_of_adduct_isotope_all_result.csv`, and
#' `clustering_of_adduct_isotope_main_result.csv`.
#'
#' @param met_raw_wide A feature table in wide format. The first column must be
#'   `mz`, the second column must be `time`, and the remaining columns are sample
#'   intensities.
#' @param metabolite_database Optional metabolite database with the columns
#'   `Mono_mass`, `Most_abundant_isotopologue_mass`, `Exact_mass`,
#'   `Exact_mass_most_abundant_isotopologue`, `Name`, `Formula`, `HMDB_ID`,
#'   `KEGG_ID`, `InChIKey`, `Subclass`, `has_NH3`, `has_H2O`, `n_NH3`, and
#'   `n_H2O`. When `NULL`, the database specified by `database` is loaded
#'   automatically.
#' @param database Database used when `metabolite_database` is `NULL`. Choose
#'   either `"metorigindb"` or `"pubchem"`. Default is `"metorigindb"`.
#' @param hmdb_concentration Optional HMDB concentration table. If `NULL`, the
#'   packaged HMDB concentration table is used via `hmdb_concentration_loading()`.
#' @param mz_threshold The ppm threshold used for m/z matching. Default is `10`.
#' @param biospecimen Biospecimen used to prioritize HMDB concentrations.
#'   Default is `"Blood"`.
#' @param All_Adduct Adducts considered during m/z matching. When `NULL`, a
#'   positive- or negative-mode default is selected from `ion_mode`.
#' @param ion_mode Ion mode string used in outputs. Default is `"positive"`.
#' @param adduct_correlation_r_threshold Spearman correlation threshold for
#'   grouping adducts. Default is `0.39`.
#' @param adduct_correlation_time_threshold Retention-time window in seconds for
#'   adduct correlation grouping. Default is `6`.
#' @param isotopic_correlation_r_threshold Spearman correlation threshold for
#'   isotope pairing. Default is `0.71`.
#' @param isotopic_correlation_time_threshold Retention-time window in seconds
#'   for isotope pairing. Default is `4`.
#' @param imputation_method Missing-value imputation method. One of
#'   `"half_min"`, `"QRILC"`, or `NA`. Default is `"half_min"`.
#' @param output_dir Directory where CSV outputs are written. Default is `"."`.
#' @param write_output Logical; if `TRUE`, write the four CSV outputs. Default
#'   is `TRUE`.
#' @param show_progress Logical; if `TRUE`, emit console progress messages and
#'   a text progress bar. Default is `TRUE`.
#' @param progress_callback Optional function for external progress updates.
#'   When provided, it is called as `progress_callback(step, total, detail)`.
#'
#' @return A named list containing the four result tables:
#'   `mz_only`, `mz_only_isotope`, `clustering_all`, and `clustering_main`.
#' @export
mz_match_clustering <- function(
    met_raw_wide,
    metabolite_database = NULL,
    database = c("metorigindb", "pubchem"),
    hmdb_concentration = NULL,
    mz_threshold = 10,
    biospecimen = "Blood",
    All_Adduct = NULL,
    ion_mode = "positive",
    adduct_correlation_r_threshold = 0.39,
    adduct_correlation_time_threshold = 6,
    isotopic_correlation_r_threshold = 0.71,
    isotopic_correlation_time_threshold = 4,
    imputation_method = "half_min",
    output_dir = ".",
    write_output = TRUE,
    show_progress = TRUE,
    progress_callback = NULL
) {
    library(dplyr)
    library(MetaboCoreUtilsAdduct)

    database <- match.arg(database)
    feature_intensity_columns <- if (ncol(met_raw_wide) > 2) {
        colnames(met_raw_wide)[-c(1, 2)]
    } else {
        character(0)
    }

    total_progress_steps <- 7L
    current_progress_step <- 0
    progress_bar <- NULL

    if (isTRUE(show_progress)) {
        progress_bar <- utils::txtProgressBar(min = 0, max = total_progress_steps, style = 3)
        on.exit(close(progress_bar), add = TRUE)
    }

    set_progress <- function(step_value, detail, announce = TRUE) {
        current_progress_step <<- max(
            current_progress_step,
            min(total_progress_steps, as.numeric(step_value))
        )
        if (isTRUE(show_progress)) {
            if (isTRUE(announce)) {
                message("[", round(current_progress_step, 2), "/", total_progress_steps, "] ", detail)
            }
            utils::setTxtProgressBar(progress_bar, current_progress_step)
        }
        if (is.function(progress_callback)) {
            try(progress_callback(current_progress_step, total_progress_steps, detail), silent = TRUE)
        }
    }

    report_progress <- function(detail) {
        set_progress(current_progress_step + 1, detail, announce = TRUE)
    }

    report_clustering_progress <- function(group_index, group_total) {
        if (group_total <= 0) {
            set_progress(total_progress_steps, "Clustering completed.", announce = TRUE)
            return(invisible(NULL))
        }

        progress_value <- 6 + (group_index / group_total)
        detail <- paste0("Clustering group ", group_index, "/", group_total)
        announce <- isTRUE(show_progress) && (
            group_index == 1 ||
                group_index == group_total ||
                group_index %% max(1, floor(group_total / 20)) == 0
        )
        set_progress(progress_value, detail, announce = announce)
    }

    default_adducts <- function(current_ion_mode) {
        if (identical(current_ion_mode, "negative")) {
            return(c("M-H", "M+Cl", "M-H-H2O", "2M-H"))
        }

        c(
            "M+H", "M+Na", "M+2Na-H", "M+H-H2O", "M+H-NH3",
            "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O"
        )
    }

    estimate_neutral_mass_range <- function(mz_values, adducts, ppm_threshold) {
        adduct_definition <- adduct_definition_loading() |>
            dplyr::filter(.data$name %in% adducts) |>
            dplyr::select(name, mass_multi, mass_add)

        if (nrow(adduct_definition) == 0) {
            stop("None of the selected adducts are available in adduct_definition.", call. = FALSE)
        }

        ppm_margin <- max(ppm_threshold, 0) / 1000000
        mz_lower <- min(mz_values) * (1 - ppm_margin)
        mz_upper <- max(mz_values) * (1 + ppm_margin)

        mass_at_mz_lower <- (mz_lower - adduct_definition$mass_add) / adduct_definition$mass_multi
        mass_at_mz_upper <- (mz_upper - adduct_definition$mass_add) / adduct_definition$mass_multi
        lower_candidates <- pmin(mass_at_mz_lower, mass_at_mz_upper)
        upper_candidates <- pmax(mass_at_mz_lower, mass_at_mz_upper)

        lower <- min(lower_candidates, na.rm = TRUE) - 1.01
        upper <- max(upper_candidates, na.rm = TRUE) + 1.01

        list(lower = lower, upper = upper)
    }

    compound_identifier_columns <- c(
        "unique_identifier", "InChIKey", "SMILES"
    )
    database_identifier_columns <- c(
        "KEGG_ID", "HMDB_ID", "CAS_ID", "INCHI_ID", "KEGG_DRUG_ID",
        "CHEMSPIDER_ID", "DRUGBANK_ID", "FOODB_ID", "PUBCHEM_COMPOUND_ID",
        "PUBCHEM_SUBSTANCE_ID", "CHEBI_ID", "CHEMBL_ID", "PDB_CCD_ID",
        "3DMET_ID", "NIKKAJI_ID", "KNAPSACK_ID", "LIPIDMAPS_ID",
        "LIPIDBANK_ID", "BIOCYC_ID", "BIGG_ID", "BIGG_IDENTIFIER_ID",
        "WIKIPEDIA_ID", "METLIN_ID", "T3DB_ID", "REACTOME_ID",
        "MODELSEED_ID", "MIMEDB_ID", "LOTUS_ID"
    )
    exact_mass_columns <- c("Exact_mass", "Exact_mass_most_abundant_isotopologue")
    annotation_detail_columns <- c(
        "Name", "Formula", "Mono_mass", "Most_abundant_isotopologue_mass",
        "Charge_natural", "Subclass",
        compound_identifier_columns,
        database_identifier_columns
    )

    validate_columns <- function(data, required_columns, label) {
        missing_columns <- setdiff(required_columns, colnames(data))
        if (length(missing_columns) > 0) {
            stop(
                label,
                " is missing required columns: ",
                paste(missing_columns, collapse = ", "),
                call. = FALSE
            )
        }
    }

    ensure_optional_columns <- function(data, character_columns = character(), numeric_columns = character()) {
        data <- dplyr::as_tibble(data)

        for (col_name in setdiff(character_columns, colnames(data))) {
            data[[col_name]] <- NA_character_
        }
        for (col_name in setdiff(numeric_columns, colnames(data))) {
            data[[col_name]] <- NA_real_
        }

        data
    }

    organize_match_output_columns <- function(data) {
        data |>
            dplyr::select(
                dplyr::any_of(c(
                    "ion_mode", "mz", "time", "mz_time", "theoretical_mz",
                    "mz_matching_ppm", "Adduct", "isotope", "mean_intensity"
                )),
                dplyr::any_of(annotation_detail_columns),
                dplyr::everything(),
                -dplyr::any_of(exact_mass_columns)
            ) |>
            dplyr::relocate(dplyr::any_of("theoretical_mz"), .before = dplyr::any_of("mz_matching_ppm"))
    }

    normalize_inchikey <- function(x) {
        toupper(trimws(as.character(x)))
    }

    clean_name <- function(x) {
        vapply(strsplit(as.character(x), ";", fixed = TRUE), `[`, character(1), 1)
    }

    prepare_hmdb_concentration_simple <- function(hmdb_concentration_input, metabolite_database_input, current_biospecimen) {
        concentration_df <- hmdb_concentration_input
        if (is.null(concentration_df)) {
            concentration_df <- hmdb_concentration_loading()
        }

        concentration_df <- dplyr::as_tibble(concentration_df)

        hmdb_id_column <- if ("HMDBID" %in% colnames(concentration_df)) {
            "HMDBID"
        } else if ("HMDB_ID" %in% colnames(concentration_df)) {
            "HMDB_ID"
        } else {
            stop(
                "hmdb_concentration must include either HMDBID or HMDB_ID.",
                call. = FALSE
            )
        }

        required_concentration_columns <- c(hmdb_id_column, "Concentration_average", "Concentration_units")
        validate_columns(concentration_df, required_concentration_columns, "hmdb_concentration")

        if ("Biospecimen" %in% colnames(concentration_df)) {
            concentration_df <- concentration_df |>
                dplyr::filter(.data$Biospecimen == current_biospecimen)
        }

        concentration_df <- concentration_df |>
            dplyr::select(
                HMDBID = dplyr::all_of(hmdb_id_column),
                Concentration_average,
                Concentration_units
            ) |>
            dplyr::filter(!is.na(.data$HMDBID)) |>
            dplyr::group_by(.data$HMDBID) |>
            dplyr::slice_max(.data$Concentration_average, n = 1, with_ties = FALSE) |>
            dplyr::ungroup()

        hmdb_mapping <- metabolite_database_input |>
            dplyr::transmute(
                HMDBID = .data$HMDB_ID,
                InChIKey = normalize_inchikey(.data$InChIKey)
            ) |>
            dplyr::filter(!is.na(.data$HMDBID), !is.na(.data$InChIKey)) |>
            dplyr::distinct()

        concentration_df |>
            dplyr::left_join(hmdb_mapping, by = "HMDBID") |>
            dplyr::filter(!is.na(.data$InChIKey)) |>
            dplyr::group_by(.data$InChIKey) |>
            dplyr::slice_max(.data$Concentration_average, n = 1, with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::select(InChIKey, Concentration_average, Concentration_units)
    }

    impute_half_min <- function(values) {
        non_missing <- values[!is.na(values)]
        if (length(non_missing) == 0) {
            return(values)
        }

        values[is.na(values)] <- min(non_missing) / 2
        values
    }

    preprocess_feature_table <- function(feature_table, method) {
        feature_table <- dplyr::as_tibble(feature_table)
        colnames(feature_table)[1:2] <- c("mz", "time")
        feature_table[feature_table == 0] <- NA

        feature_table <- feature_table |>
            dplyr::mutate(mean_intensity = rowMeans(dplyr::pick(-c("mz", "time")), na.rm = TRUE)) |>
            dplyr::group_by(.data$mz, .data$time) |>
            dplyr::slice_max(.data$mean_intensity, n = 1, with_ties = FALSE) |>
            dplyr::ungroup()

        feature_table_original_mean_intensity <- feature_table |>
            dplyr::select(mz, time, mean_intensity)

        feature_table <- feature_table |>
            dplyr::select(-mean_intensity)

        transformed <- feature_table
        transformed[, -c(1, 2)] <- log2(as.data.frame(transformed[, -c(1, 2)]))
        transformed_matrix <- as.matrix(transformed[, -c(1, 2)])
        transformed_matrix[is.nan(transformed_matrix)] <- NA_real_
        transformed[, -c(1, 2)] <- transformed_matrix

        if (identical(method, "QRILC")) {
            meta <- transformed[, 1:2]
            intensity <- as.data.frame(transformed[, -c(1, 2)])
            intensity_imputed <- imputeLCMD::impute.QRILC(intensity)[[1]]
            transformed <- dplyr::bind_cols(meta, dplyr::as_tibble(intensity_imputed))
            colnames(transformed) <- c("mz", "time", colnames(intensity))
        } else if (identical(method, "half_min")) {
            meta <- transformed[, 1:2]
            intensity <- as.data.frame(transformed[, -c(1, 2)])
            intensity_imputed <- t(apply(intensity, 1, impute_half_min))
            intensity_imputed <- dplyr::as_tibble(as.data.frame(intensity_imputed))
            colnames(intensity_imputed) <- colnames(intensity)
            transformed <- dplyr::bind_cols(meta, intensity_imputed)
        } else if (!(length(method) == 1 && is.na(method))) {
            stop("imputation_method must be 'half_min', 'QRILC', or NA.", call. = FALSE)
        }

        normalized_intensity <- preprocessCore::normalize.quantiles(
            as.matrix(dplyr::select(transformed, -c(mz, time)))
        )
        colnames(normalized_intensity) <- colnames(transformed)[-c(1, 2)]
        normalized <- dplyr::bind_cols(
            dplyr::select(transformed, mz, time),
            dplyr::as_tibble(as.data.frame(normalized_intensity))
        )

        list(
            processed = normalized,
            mean_intensity = feature_table_original_mean_intensity
        )
    }

    write_result_csv <- function(data, file_name) {
        if (isTRUE(write_output)) {
            readr::write_csv(data, file.path(output_dir, file_name))
        }
    }

    drop_feature_intensity_columns <- function(data, intensity_cols) {
        dplyr::as_tibble(data) |>
            dplyr::select(-dplyr::any_of(intensity_cols))
    }

    if (is.null(All_Adduct)) {
        All_Adduct <- default_adducts(ion_mode)
    }

    if (!dir.exists(output_dir) && isTRUE(write_output)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    met_raw_wide <- dplyr::as_tibble(met_raw_wide)
    if (ncol(met_raw_wide) < 3) {
        stop("met_raw_wide must contain mz, time, and at least one sample column.", call. = FALSE)
    }
    colnames(met_raw_wide)[1:2] <- c("mz", "time")

    mz_values <- suppressWarnings(as.numeric(met_raw_wide$mz))
    mz_values <- mz_values[is.finite(mz_values)]
    if (length(mz_values) == 0) {
        stop("No finite mz values found in met_raw_wide.", call. = FALSE)
    }

    if (is.null(metabolite_database)) {
        report_progress(paste0("Loading ", database, " database..."))
    } else {
        report_progress("Using provided metabolite database...")
    }

    if (is.null(metabolite_database)) {
        neutral_mass_range <- estimate_neutral_mass_range(
            mz_values = mz_values,
            adducts = All_Adduct,
            ppm_threshold = mz_threshold
        )

        metabolite_database <- metabolite_database_loading(database = database) |>
            dplyr::filter(
                (.data$Exact_mass >= neutral_mass_range$lower & .data$Exact_mass <= neutral_mass_range$upper) |
                    (
                        .data$Exact_mass_most_abundant_isotopologue >= neutral_mass_range$lower &
                            .data$Exact_mass_most_abundant_isotopologue <= neutral_mass_range$upper
                    )
            ) |>
            dplyr::collect()
    }

    required_metabolite_columns <- c(
        "Mono_mass", "Most_abundant_isotopologue_mass",
        "Exact_mass", "Exact_mass_most_abundant_isotopologue",
        "Name", "Formula",
        "HMDB_ID", "KEGG_ID", "InChIKey", "Subclass", "has_NH3", "has_H2O",
        "n_NH3", "n_H2O"
    )
    validate_columns(metabolite_database, required_metabolite_columns, "metabolite_database")

    metabolite_database <- dplyr::as_tibble(metabolite_database) |>
        dplyr::mutate(InChIKey = normalize_inchikey(.data$InChIKey)) |>
        dplyr::group_by(.data$InChIKey) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
    metabolite_database <- ensure_optional_columns(
        metabolite_database,
        character_columns = unique(c("Formula", compound_identifier_columns, database_identifier_columns)),
        numeric_columns = c("Charge_natural")
    )

    report_progress("Preparing HMDB concentration mapping...")
    hmdb_metabolites_concentrations_average_simple <- prepare_hmdb_concentration_simple(
        hmdb_concentration_input = hmdb_concentration,
        metabolite_database_input = metabolite_database,
        current_biospecimen = biospecimen
    )

    report_progress("Building adduct/isotope m/z library...")
    metabolite_mz_library <- tidyr::crossing(
        metabolite_database,
        Adduct = All_Adduct
    ) |>
        dplyr::filter(
            (.data$Adduct == "M+H-NH3" & .data$has_NH3 == TRUE) |
                (.data$Adduct %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & .data$n_H2O >= 1) |
                (.data$Adduct %in% c("M-2H2O+H", "M+H-2H2O") & .data$n_H2O >= 2) |
                (!.data$Adduct %in% c(
                    "M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H",
                    "M-H-H2O", "M-2H2O+H", "M+H-2H2O"
                ))
        )

    metabolite_mz_library$mz <- mass2mz_df_safe(
        mass = metabolite_mz_library$Exact_mass,
        adduct = metabolite_mz_library$Adduct
    )$mz
    metabolite_mz_library$mz_isotope <- mass2mz_df_safe(
        mass = metabolite_mz_library$Exact_mass_most_abundant_isotopologue,
        adduct = metabolite_mz_library$Adduct
    )$mz
    metabolite_mz_library <- metabolite_mz_library |>
        dplyr::filter(.data$mz > 0)

    report_progress("Preprocessing feature table...")
    preprocessing <- preprocess_feature_table(met_raw_wide, imputation_method)
    met_raw_wide_processed <- preprocessing$processed
    met_raw_wide_original_mean_intensity <- preprocessing$mean_intensity

    mz_matching <- function(feature_table, mz_library) {
        overlap <- find.Overlapping.mzs(
            data.frame(mz = feature_table$mz, stringsAsFactors = FALSE),
            data.frame(mz = mz_library$mz, stringsAsFactors = FALSE),
            mz.thresh = mz_threshold
        )

        if (nrow(overlap) == 0) {
            return(list(internal = dplyr::tibble(), output = dplyr::tibble()))
        }

        metabolite_matched <- dplyr::slice(mz_library, overlap$index.B) |>
            dplyr::rename(mz_annotated = mz)
        feature_matched <- dplyr::slice(feature_table, overlap$index.A)
        colnames(feature_matched)[1:2] <- c("mz_sample", "time_sample")

        annotated <- dplyr::bind_cols(metabolite_matched, feature_matched) |>
            dplyr::left_join(
                met_raw_wide_original_mean_intensity,
                by = c("mz_sample" = "mz", "time_sample" = "time")
            ) |>
            dplyr::relocate(mean_intensity, .after = time_sample) |>
            dplyr::mutate(
                mean_intensity = round(.data$mean_intensity),
                ion_mode = ion_mode,
                mz_sample = round(.data$mz_sample, 4),
                time_sample = round(.data$time_sample, 0),
                mz_time_sample = paste0(.data$mz_sample, "_", .data$time_sample),
                isotope = NA_character_,
                Mono_mass = round(.data$Mono_mass, 4)
            ) |>
            dplyr::select(
                ion_mode, KEGG_ID, HMDB_ID, mz_annotated, Name,
                Formula, Mono_mass, Adduct, isotope,
                mz_sample, time_sample, mz_time_sample,
                mean_intensity, dplyr::everything(), -c(Most_abundant_isotopologue_mass, mz_isotope)
            ) |>
            dplyr::select(-dplyr::any_of(c("Mass_diff", "CID"))) |>
            dplyr::distinct()

        # Safety guard: retain only matches within the requested ppm threshold.
        # This protects downstream results even if overlap backends return spurious pairs.
        annotated <- annotated |>
            dplyr::mutate(
                mz_matching_ppm = abs((.data$mz_sample - .data$mz_annotated) / .data$mz_annotated * 1000000)
            ) |>
            dplyr::filter(.data$mz_matching_ppm <= mz_threshold + 1e-12)

        if (nrow(annotated) == 0) {
            return(list(internal = dplyr::tibble(), output = dplyr::tibble()))
        }

        output <- annotated |>
            dplyr::rename(
                theoretical_mz = mz_annotated,
                mz = mz_sample,
                time = time_sample,
                mz_time = mz_time_sample
            ) |>
            dplyr::mutate(
                theoretical_mz = round(.data$theoretical_mz, 4),
                mz = round(.data$mz, 4),
                time = round(.data$time, 0),
                mz_matching_ppm = round(.data$mz_matching_ppm, 1),
                Name = clean_name(.data$Name)
            ) |>
            dplyr::relocate(theoretical_mz, .before = mz_matching_ppm) |>
            dplyr::relocate(mz_matching_ppm, .after = mz_time) |>
            organize_match_output_columns()

        list(internal = annotated, output = output)
    }

    mz_isotope_matching <- function(feature_table, mz_library) {
        overlap <- find.Overlapping.mzs(
            data.frame(mz = feature_table$mz, stringsAsFactors = FALSE),
            data.frame(mz = mz_library$mz_isotope, stringsAsFactors = FALSE),
            mz.thresh = mz_threshold
        )

        if (nrow(overlap) == 0) {
            return(list(internal = dplyr::tibble(), output = dplyr::tibble()))
        }

        metabolite_matched <- dplyr::slice(mz_library, overlap$index.B) |>
            dplyr::rename(mz_annotated_isotope = mz_isotope) |>
            dplyr::select(-mz)
        if (!"Mass_diff" %in% colnames(metabolite_matched)) {
            metabolite_matched$Mass_diff <- NA_character_
        }
        feature_matched <- dplyr::slice(feature_table, overlap$index.A)
        colnames(feature_matched)[1:2] <- c("mz_sample_isotope", "time_sample_isotope")

        annotated <- dplyr::bind_cols(metabolite_matched, feature_matched) |>
            dplyr::left_join(
                met_raw_wide_original_mean_intensity,
                by = c("mz_sample_isotope" = "mz", "time_sample_isotope" = "time")
            ) |>
            dplyr::relocate(mean_intensity, .after = time_sample_isotope) |>
            dplyr::mutate(
                mean_intensity = round(.data$mean_intensity),
                ion_mode = ion_mode,
                mz_sample_isotope = round(.data$mz_sample_isotope, 4),
                time_sample_isotope = round(.data$time_sample_isotope, 0),
                mz_time_sample_isotope = paste0(.data$mz_sample_isotope, "_", .data$time_sample_isotope),
                isotope = paste0("[", .data$Mass_diff, "]"),
                Mono_mass = round(.data$Mono_mass, 4)
            ) |>
            dplyr::select(
                ion_mode, KEGG_ID, HMDB_ID, mz_annotated_isotope,
                Name, Formula, Most_abundant_isotopologue_mass,
                Mono_mass, Adduct, isotope, mz_sample_isotope,
                time_sample_isotope, mz_time_sample_isotope,
                mean_intensity, dplyr::everything()
            ) |>
            dplyr::select(-dplyr::any_of(c("Mass_diff", "CID"))) |>
            dplyr::distinct()

        annotated <- annotated |>
            dplyr::mutate(
                mz_matching_ppm = abs((.data$mz_sample_isotope - .data$mz_annotated_isotope) / .data$mz_annotated_isotope * 1000000)
            ) |>
            dplyr::filter(.data$mz_matching_ppm <= mz_threshold + 1e-12)

        if (nrow(annotated) == 0) {
            return(list(internal = dplyr::tibble(), output = dplyr::tibble()))
        }

        output <- annotated |>
            dplyr::rename(
                theoretical_mz = mz_annotated_isotope,
                mz = mz_sample_isotope,
                time = time_sample_isotope,
                mz_time = mz_time_sample_isotope
            ) |>
            dplyr::mutate(
                theoretical_mz = round(.data$theoretical_mz, 4),
                mz = round(.data$mz, 4),
                time = round(.data$time, 0),
                Most_abundant_isotopologue_mass = round(.data$Most_abundant_isotopologue_mass, 4),
                mz_matching_ppm = round(.data$mz_matching_ppm, 1),
                Name = clean_name(.data$Name)
            ) |>
            dplyr::relocate(theoretical_mz, .before = mz_matching_ppm) |>
            dplyr::relocate(mz_matching_ppm, .after = mz_time) |>
            organize_match_output_columns()

        list(internal = annotated, output = output)
    }

    report_progress("Running m/z matching...")
    mz_only_result <- mz_matching(met_raw_wide_processed, metabolite_mz_library)
    report_progress("Running m/z isotopologue matching...")
    mz_only_isotope_result <- mz_isotope_matching(met_raw_wide_processed, metabolite_mz_library)

    mz_only_internal <- mz_only_result$internal
    mz_only_output <- mz_only_result$output
    mz_only_isotope_internal <- mz_only_isotope_result$internal
    mz_only_isotope_output <- mz_only_isotope_result$output

    mz_only_output <- drop_feature_intensity_columns(mz_only_output, feature_intensity_columns)
    mz_only_isotope_output <- drop_feature_intensity_columns(mz_only_isotope_output, feature_intensity_columns)

    write_result_csv(
        mz_only_output,
        paste0(ion_mode, "_metabolite_annotated_", mz_threshold, "ppm_mzonly.csv")
    )
    write_result_csv(
        mz_only_isotope_output,
        paste0(ion_mode, "_metabolite_annotated_", mz_threshold, "ppm_mzonly_isotope.csv")
    )

    if (nrow(mz_only_internal) == 0) {
        set_progress(
            total_progress_steps,
            "No primary m/z matches found for the selected settings. Skipping clustering.",
            announce = TRUE
        )

        clustering_all_empty <- dplyr::tibble(
            Mono_mass = numeric(),
            identification_type_annotated = character(),
            ion_mode = character(),
            Adduct_annotated = character(),
            mz_annotated = numeric(),
            time_annotated = numeric(),
            mz_time_annotated = character(),
            adduct_corr_cluster = logical(),
            mz_isotope = numeric(),
            time_isotope = numeric(),
            mz_time_isotope_annotated = character(),
            mean_intensity = numeric(),
            mean_intensity_isotope = numeric(),
            isotopic_correlation = numeric(),
            Probability = numeric(),
            cluster_identification = integer(),
            Name = character(),
            Formula = character(),
            Most_abundant_isotopologue_mass = numeric(),
            Charge_natural = numeric(),
            Subclass = character(),
            has_NH3 = logical(),
            has_H2O = logical(),
            n_NH3 = integer(),
            n_H2O = integer(),
            Concentration_average = numeric(),
            Concentration_units = character()
        )
        clustering_all_empty <- ensure_optional_columns(
            clustering_all_empty,
            character_columns = c(compound_identifier_columns, database_identifier_columns)
        ) |>
            dplyr::select(
                dplyr::any_of(c(
                    "Mono_mass", "identification_type_annotated", "ion_mode",
                    "Adduct_annotated", "mz_annotated", "time_annotated",
                    "mz_time_annotated", "mean_intensity", "adduct_corr_cluster",
                    "mz_isotope", "time_isotope", "mz_time_isotope_annotated",
                    "mean_intensity_isotope", "isotopic_correlation",
                    "Probability", "cluster_identification"
                )),
                dplyr::any_of(annotation_detail_columns),
                dplyr::any_of(c(
                    "Concentration_average", "Concentration_units",
                    "has_NH3", "has_H2O", "n_NH3", "n_H2O"
                )),
                dplyr::everything()
            )

        clustering_main_empty <- dplyr::tibble(
            Mono_mass = numeric(),
            identification_type = character(),
            ion_mode = character(),
            identified_Name = character(),
            Formula = character(),
            Most_abundant_isotopologue_mass = numeric(),
            Charge_natural = numeric(),
            Subclass = character(),
            Adduct = character(),
            mz = numeric(),
            time = numeric(),
            mz_time = character(),
            cluster_identification = integer(),
            identification_method = character(),
            Concentration_average = numeric(),
            Concentration_units = character(),
            Probability = numeric(),
            match_category = character()
        )
        clustering_main_empty <- ensure_optional_columns(
            clustering_main_empty,
            character_columns = c(compound_identifier_columns, database_identifier_columns)
        ) |>
            dplyr::select(
                dplyr::any_of(c(
                    "Mono_mass", "identification_type", "ion_mode", "identified_Name",
                    "Formula", "Most_abundant_isotopologue_mass",
                    "Charge_natural", "Subclass"
                )),
                dplyr::any_of(c(compound_identifier_columns, database_identifier_columns)),
                dplyr::any_of(c(
                    "Adduct", "mz", "time", "mz_time", "cluster_identification",
                    "identification_method", "Concentration_average",
                    "Concentration_units", "Probability", "match_category"
                )),
                dplyr::everything()
            )

        write_result_csv(clustering_all_empty, "clustering_of_adduct_isotope_all_result.csv")
        write_result_csv(clustering_main_empty, "clustering_of_adduct_isotope_main_result.csv")

        return(list(
            mz_only = mz_only_output,
            mz_only_isotope = mz_only_isotope_output,
            clustering_all = clustering_all_empty,
            clustering_main = clustering_main_empty
        ))
    }

    met_raw_wide_processed$mz <- round(met_raw_wide_processed$mz, 4)
    met_raw_wide_processed$time <- round(met_raw_wide_processed$time, 0)
    met_raw_long <- met_raw_wide_processed |>
        dplyr::mutate(mz_time = paste0(.data$mz, "_", .data$time)) |>
        dplyr::relocate(mz_time, .after = time) |>
        dplyr::mutate(mean_intensity = rowMeans(dplyr::pick(-c("mz", "time", "mz_time")), na.rm = TRUE)) |>
        dplyr::group_by(.data$mz_time) |>
        dplyr::slice_max(.data$mean_intensity, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-mean_intensity)

    correlation_input <- as.data.frame(t(met_raw_long[, 4:ncol(met_raw_long)]))
    colnames(correlation_input) <- met_raw_long$mz_time
    correlation_input <- dplyr::as_tibble(correlation_input)

    met_raw_wide_final_monomass <- mz_only_internal |>
        dplyr::mutate(
            InChIKey_mz_time = paste0(.data$InChIKey, "_", .data$mz_sample, "_", .data$time_sample)
        ) |>
        dplyr::group_by(.data$InChIKey_mz_time) |>
        dplyr::slice_max(.data$mean_intensity, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-InChIKey_mz_time) |>
        dplyr::mutate(identification_type_annotated = "identified") |>
        dplyr::rename(Confirmed_Name = Name) |>
        dplyr::select(
            Exact_mass, Mono_mass, InChIKey, identification_type_annotated,
            ion_mode, Confirmed_Name, Adduct, mz_sample,
            time_sample, mz_time_sample
        ) |>
        dplyr::rename(
            Adduct_annotated = Adduct,
            mz_annotated = mz_sample,
            time_annotated = time_sample,
            mz_time_annotated = mz_time_sample
        ) |>
        dplyr::mutate(Mono_mass = round(.data$Mono_mass, 4))

    met_raw_wide_final_monomass_isotope <- mz_only_isotope_internal |>
        dplyr::mutate(
            InChIKey_mz_time = paste0(.data$InChIKey, "_", .data$mz_sample_isotope, "_", .data$time_sample_isotope)
        ) |>
        dplyr::group_by(.data$InChIKey_mz_time) |>
        dplyr::slice_max(.data$mean_intensity, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-InChIKey_mz_time) |>
        dplyr::mutate(identification_type_annotated = "identified") |>
        dplyr::rename(Confirmed_Name = Name) |>
        dplyr::select(
            Exact_mass, Mono_mass, InChIKey, identification_type_annotated,
            ion_mode, Confirmed_Name, Adduct, mz_sample_isotope,
            time_sample_isotope, mz_time_sample_isotope
        ) |>
        dplyr::rename(
            Adduct_annotated = Adduct,
            mz_annotated = mz_sample_isotope,
            time_annotated = time_sample_isotope,
            mz_time_annotated = mz_time_sample_isotope
        ) |>
        dplyr::mutate(Mono_mass = round(.data$Mono_mass, 4))

    met_raw_wide_final_monomass_simple <- met_raw_wide_final_monomass |>
        dplyr::arrange(.data$Adduct_annotated, .data$Exact_mass) |>
        dplyr::select(-InChIKey, -Confirmed_Name) |>
        dplyr::distinct()

    met_raw_wide_final_monomass_isotope_simple <- met_raw_wide_final_monomass_isotope |>
        dplyr::arrange(.data$Adduct_annotated, .data$Exact_mass) |>
        dplyr::select(-InChIKey, -Confirmed_Name) |>
        dplyr::distinct()

    if (isTRUE(show_progress)) {
        message("[6/7] Running clustering...")
    }
    if (is.function(progress_callback)) {
        try(progress_callback(current_progress_step, total_progress_steps, "Running clustering..."), silent = TRUE)
    }
    data_split_cluster <- split(met_raw_wide_final_monomass_simple, met_raw_wide_final_monomass_simple$Exact_mass)

    bayesian_probability_calculation_cluster <- function(data) {
        current_exact_mass <- unique(data$Exact_mass)

        find_high_correlation_pairs <- function(cor_matrix, threshold) {
            cor_matrix[lower.tri(cor_matrix)] <- NA
            diag(cor_matrix) <- NA
            high_cor_indices <- which(cor_matrix >= threshold, arr.ind = TRUE)

            if (length(high_cor_indices) == 0) {
                return(list())
            }

            high_cor_pairs <- apply(high_cor_indices, 1, function(idx) {
                colnames(cor_matrix)[idx]
            })

            unique(apply(high_cor_pairs, 2, function(pair) {
                sort(pair)
            }, simplify = FALSE))
        }

        if (nrow(data) > 1) {
            data_ordered <- data |>
                dplyr::arrange(.data$time_annotated)

            time_groups <- list()
            current_group <- data_ordered$mz_time_annotated[1]
            last_time <- data_ordered$time_annotated[1]

            if (nrow(data_ordered) > 1) {
                for (i in 2:nrow(data_ordered)) {
                    if ((data_ordered$time_annotated[i] - last_time) <= adduct_correlation_time_threshold) {
                        current_group <- c(current_group, data_ordered$mz_time_annotated[i])
                    } else {
                        if (length(current_group) > 1) {
                            time_groups[[length(time_groups) + 1]] <- current_group
                        }
                        current_group <- data_ordered$mz_time_annotated[i]
                    }
                    last_time <- data_ordered$time_annotated[i]
                }
            }
            if (length(current_group) > 1) {
                time_groups[[length(time_groups) + 1]] <- current_group
            }

            filtered_original_groups <- lapply(time_groups, function(group) {
                cor_input <- correlation_input[, group, drop = FALSE]
                cor_matrix <- suppressWarnings(stats::cor(
                    cor_input,
                    cor_input,
                    method = "spearman",
                    use = "pairwise.complete.obs"
                ))
                find_high_correlation_pairs(cor_matrix, adduct_correlation_r_threshold)
            })

            flat_list <- do.call(c, filtered_original_groups)

            if (length(flat_list) > 0) {
                intra_correlation_all_clean_df <- data.frame(
                    adduct_corr_cluster = rep(seq_along(flat_list), lengths(flat_list)),
                    variable = unlist(flat_list)
                )

                data_adduct_corr <- data |>
                    dplyr::full_join(
                        intra_correlation_all_clean_df,
                        by = c("mz_time_annotated" = "variable")
                    )
            } else {
                data_adduct_corr <- data |>
                    dplyr::mutate(adduct_corr_cluster = NA)
            }
        } else {
            data_adduct_corr <- data |>
                dplyr::mutate(adduct_corr_cluster = NA)
        }

        data_adduct_corr_2 <- data_adduct_corr |>
            dplyr::distinct(.data$Exact_mass, .data$Adduct_annotated, .data$mz_annotated, .data$time_annotated, .keep_all = TRUE) |>
            dplyr::mutate(adduct_corr_cluster = ifelse(!is.na(.data$adduct_corr_cluster), TRUE, NA))

        abundance_ratios <- data.frame()
        primary_df <- met_raw_wide_final_monomass_simple |>
            dplyr::filter(.data$Exact_mass %in% current_exact_mass)
        isotope_df <- met_raw_wide_final_monomass_isotope_simple |>
            dplyr::filter(.data$Exact_mass %in% current_exact_mass)

        if (nrow(isotope_df) == 0) {
            data_adduct_corr_3 <- data_adduct_corr_2 |>
                dplyr::mutate(mz_isotope = NA_real_, time_isotope = NA_real_)
        } else if (any(primary_df$Adduct_annotated %in% isotope_df$Adduct_annotated)) {
            for (i in seq_len(nrow(primary_df))) {
                current_primary_df <- primary_df[i, ]
                isotope_df_filtered <- isotope_df |>
                    dplyr::filter(.data$Adduct_annotated == current_primary_df$Adduct_annotated) |>
                    dplyr::filter(
                        .data$time_annotated >= current_primary_df$time_annotated - isotopic_correlation_time_threshold,
                        .data$time_annotated <= current_primary_df$time_annotated + isotopic_correlation_time_threshold
                    )

                if (nrow(isotope_df_filtered) == 0) {
                    next
                }

                for (j in seq_len(nrow(isotope_df_filtered))) {
                    current_isotope_df <- isotope_df_filtered[j, ]
                    primary_adduct_intensity <- correlation_input[, current_primary_df$mz_time_annotated, drop = TRUE]
                    isotope_adduct_intensity <- correlation_input[, current_isotope_df$mz_time_annotated, drop = TRUE]
                    ratios <- 2^(isotope_adduct_intensity - primary_adduct_intensity)

                    result_df <- data.frame(
                        current_exact_mass = current_exact_mass,
                        mz_primary = current_primary_df$mz_annotated,
                        time_primary = current_primary_df$time_annotated,
                        mz_isotope = current_isotope_df$mz_annotated,
                        time_isotope = current_isotope_df$time_annotated,
                        Adduct = current_primary_df$Adduct_annotated,
                        primary = primary_adduct_intensity,
                        isotope = isotope_adduct_intensity,
                        ratios = ratios
                    )

                    abundance_ratios <- rbind(abundance_ratios, result_df)
                }
            }

            if (nrow(abundance_ratios) == 0) {
                data_adduct_corr_3 <- data_adduct_corr_2 |>
                    dplyr::mutate(mz_isotope = NA_real_, time_isotope = NA_real_)
            } else {
                abundance_ratios_2 <- dplyr::as_tibble(abundance_ratios)
                abundance_ratios_2$ratios[!is.finite(abundance_ratios_2$ratios)] <- NA_real_
                abundance_ratios_2$ratios[abundance_ratios_2$ratios == 0] <- NA_real_
                abundance_ratios_2$ratios <- abundance_ratios_2$ratios * 100
                abundance_ratios_2$ratios[abundance_ratios_2$ratios == 0] <- NA_real_

                abundance_ratios_4 <- abundance_ratios_2 |>
                    dplyr::group_by(
                        .data$current_exact_mass, .data$mz_primary, .data$time_primary,
                        .data$mz_isotope, .data$time_isotope, .data$Adduct
                    ) |>
                    dplyr::mutate(
                        correlation = suppressWarnings(
                            stats::cor(
                                .data$primary,
                                .data$isotope,
                                method = "spearman",
                                use = "pairwise.complete.obs"
                            )
                        )
                    ) |>
                    dplyr::summarise(
                        mean_abundance_ratio = mean(.data$ratios, na.rm = TRUE),
                        mean_correlation = mean(.data$correlation, na.rm = TRUE),
                        .groups = "drop"
                    ) |>
                    dplyr::filter(
                        .data$mean_abundance_ratio <= 100,
                        .data$mean_correlation >= isotopic_correlation_r_threshold
                    )

                if (nrow(abundance_ratios_4) > 0) {
                    abundance_ratios_4 <- abundance_ratios_4 |>
                        dplyr::transmute(
                            mz_primary = round(.data$mz_primary, 4),
                            time_primary = round(.data$time_primary, 0),
                            current_exact_mass = .data$current_exact_mass,
                            mz_isotope = .data$mz_isotope,
                            time_isotope = .data$time_isotope,
                            Adduct = .data$Adduct,
                            correlation = round(.data$mean_correlation, 2)
                        )

                    data_adduct_corr_2$mz_annotated <- round(data_adduct_corr_2$mz_annotated, 4)
                    data_adduct_corr_2$time_annotated <- round(data_adduct_corr_2$time_annotated, 0)

                    data_adduct_corr_3 <- data_adduct_corr_2 |>
                        dplyr::left_join(
                            abundance_ratios_4,
                            by = c(
                                "mz_annotated" = "mz_primary",
                                "time_annotated" = "time_primary",
                                "Exact_mass" = "current_exact_mass",
                                "Adduct_annotated" = "Adduct"
                            )
                        ) |>
                        dplyr::filter(!is.na(.data$mz_annotated))
                } else {
                    data_adduct_corr_3 <- data_adduct_corr_2 |>
                        dplyr::mutate(mz_isotope = NA_real_, time_isotope = NA_real_)
                }
            }
        } else {
            data_adduct_corr_3 <- data_adduct_corr_2 |>
                dplyr::mutate(mz_isotope = NA_real_, time_isotope = NA_real_)
        }

        data_adduct_corr_3 <- data_adduct_corr_3 |>
            dplyr::distinct(.data$Exact_mass, .data$Adduct_annotated, .data$mz_annotated, .data$time_annotated, .keep_all = TRUE)

        if (all(is.na(data_adduct_corr_3$adduct_corr_cluster)) && all(is.na(data_adduct_corr_3$mz_isotope))) {
            return(data.frame())
        }

        data_adduct_corr_3 <- data_adduct_corr_3 |>
            dplyr::filter(!is.na(.data$adduct_corr_cluster) | !is.na(.data$mz_isotope))

        if (nrow(data_adduct_corr_3) > 10) {
            data_adduct_corr_3_filtered <- data_adduct_corr_3 |>
                dplyr::filter(!is.na(.data$adduct_corr_cluster) & !is.na(.data$mz_isotope))

            if (nrow(data_adduct_corr_3_filtered) == 0) {
                data_adduct_corr_3_filtered <- data_adduct_corr_3 |>
                    dplyr::arrange(dplyr::desc(.data$Adduct_annotated %in% c("M+H", "M-H"))) |>
                    dplyr::slice_head(n = 10)
            }

            data_adduct_corr_3 <- data_adduct_corr_3_filtered
        }

        data_adduct_corr_3 <- data_adduct_corr_3 |>
            dplyr::distinct(.data$Exact_mass, .data$mz_time_annotated, .keep_all = TRUE)

        if (nrow(data_adduct_corr_3) == 1) {
            return(data_adduct_corr_3 |>
                dplyr::mutate(cluster_identification = 1L, Probability = 100))
        }

        bayesian_probability_calculation_cluster_child <- function(
            data_cluster,
            n_MH_right = 1259,
            n_MH_total = 1883,
            n_other_right = 1144,
            n_other_total = 6361,
            n_cor_yes_right = 1156,
            n_cor_yes_total = 2121,
            n_cor_no_right = 1235,
            n_cor_no_total = 6099,
            n_iso_yes_right = 537,
            n_iso_yes_total = 846,
            n_iso_no_right = 2431,
            n_iso_no_total = 10149,
            eps = 1e-12
        ) {
            clamp <- function(p) {
                pmin(pmax(p, eps), 1 - eps)
            }

            k <- nrow(data_cluster)
            n_MH_wrong <- n_MH_total - n_MH_right
            n_other_wrong <- n_other_total - n_other_right
            n_cor_yes_wrong <- n_cor_yes_total - n_cor_yes_right
            n_cor_no_wrong <- n_cor_no_total - n_cor_no_right
            n_iso_yes_wrong <- n_iso_yes_total - n_iso_yes_right
            n_iso_no_wrong <- n_iso_no_total - n_iso_no_right

            p_A_MH_given_H <- clamp(n_MH_right / (n_MH_right + n_other_right))
            p_A_other_given_H <- clamp(n_other_right / (n_MH_right + n_other_right))
            p_A_MH_given_notH <- clamp(n_MH_wrong / (n_MH_wrong + n_other_wrong))
            p_A_other_given_notH <- clamp(n_other_wrong / (n_MH_wrong + n_other_wrong))

            p_C_yes_given_H <- clamp(n_cor_yes_right / (n_cor_yes_right + n_cor_no_right))
            p_C_no_given_H <- clamp(n_cor_no_right / (n_cor_yes_right + n_cor_no_right))
            p_C_yes_given_notH <- clamp(n_cor_yes_wrong / (n_cor_yes_wrong + n_cor_no_wrong))
            p_C_no_given_notH <- clamp(n_cor_no_wrong / (n_cor_yes_wrong + n_cor_no_wrong))

            p_I_yes_given_H <- clamp(n_iso_yes_right / (n_iso_yes_right + n_iso_no_right))
            p_I_no_given_H <- clamp(n_iso_no_right / (n_iso_yes_right + n_iso_no_right))
            p_I_yes_given_notH <- clamp(n_iso_yes_wrong / (n_iso_yes_wrong + n_iso_no_wrong))
            p_I_no_given_notH <- clamp(n_iso_no_wrong / (n_iso_yes_wrong + n_iso_no_wrong))

            A_is_MH <- data_cluster$Adduct_annotated %in% c("M+H", "M-H")
            C_yes <- !is.na(data_cluster$adduct_corr_cluster) & data_cluster$adduct_corr_cluster == TRUE
            I_yes <- !is.na(data_cluster$mz_isotope)
            prior <- rep(1 / k, k)

            log_A <- ifelse(
                A_is_MH,
                log(p_A_MH_given_H / p_A_MH_given_notH),
                log(p_A_other_given_H / p_A_other_given_notH)
            )
            log_C <- ifelse(
                C_yes,
                log(p_C_yes_given_H / p_C_yes_given_notH),
                log(p_C_no_given_H / p_C_no_given_notH)
            )
            log_I <- ifelse(
                I_yes,
                log(p_I_yes_given_H / p_I_yes_given_notH),
                log(p_I_no_given_H / p_I_no_given_notH)
            )

            log_w <- log(prior) + log_A + log_C + log_I
            centered <- log_w - max(log_w, na.rm = TRUE)
            weights <- exp(centered)
            posterior <- weights / sum(weights, na.rm = TRUE)

            out <- data_cluster
            out$Probability <- 100 * posterior
            max_prob <- max(out$Probability, na.rm = TRUE)
            out$cluster_identification <- as.integer(abs(out$Probability - max_prob) < 1e-8)
            out[order(-out$Probability), , drop = FALSE]
        }

        bayesian_probability_calculation_cluster_child(data_adduct_corr_3)
    }

    results_cluster <- vector("list", length(data_split_cluster))
    if (length(data_split_cluster) == 0) {
        set_progress(total_progress_steps, "Clustering completed.", announce = TRUE)
    } else {
        for (i in seq_along(data_split_cluster)) {
            results_cluster[[i]] <- bayesian_probability_calculation_cluster(data_split_cluster[[i]])
            report_clustering_progress(i, length(data_split_cluster))
        }
    }
    final_results_cluster <- dplyr::bind_rows(results_cluster)

    if (nrow(final_results_cluster) == 0 && ncol(final_results_cluster) == 0) {
        final_results_cluster <- dplyr::tibble(
            Exact_mass = numeric(),
            Mono_mass = numeric(),
            identification_type_annotated = character(),
            ion_mode = character(),
            Adduct_annotated = character(),
            mz_annotated = numeric(),
            time_annotated = numeric(),
            mz_time_annotated = character(),
            adduct_corr_cluster = logical(),
            mz_isotope = numeric(),
            time_isotope = numeric(),
            correlation = numeric(),
            Probability = numeric(),
            cluster_identification = integer()
        )
    }

    final_results_cluster <- dplyr::as_tibble(final_results_cluster) |>
        dplyr::rename(isotopic_correlation = dplyr::any_of("correlation"))

    met_raw_wide_original_mean_intensity_2 <- met_raw_wide_original_mean_intensity |>
        dplyr::mutate(
            mz_time = paste0(round(.data$mz, 4), "_", round(.data$time, 0)),
            mean_intensity = round(.data$mean_intensity, 0)
        ) |>
        dplyr::select(mz_time, mean_intensity)

    final_results_cluster <- final_results_cluster |>
        dplyr::mutate(
            mz_time_isotope_annotated = paste0(round(.data$mz_isotope, 4), "_", round(.data$time_isotope, 0))
        ) |>
        dplyr::relocate(mz_time_isotope_annotated, .after = time_isotope) |>
        dplyr::left_join(
            met_raw_wide_original_mean_intensity_2,
            by = c("mz_time_annotated" = "mz_time")
        ) |>
        dplyr::left_join(
            met_raw_wide_original_mean_intensity_2,
            by = c("mz_time_isotope_annotated" = "mz_time")
        ) |>
        dplyr::rename(
            mean_intensity = mean_intensity.x,
            mean_intensity_isotope = mean_intensity.y
        ) |>
        dplyr::relocate(mean_intensity, .after = mz_time_annotated) |>
        dplyr::relocate(mean_intensity_isotope, .after = mz_time_isotope_annotated)

    # remove self-correlation rows (same feature annotated as both primary and isotope)
    final_results_cluster <- final_results_cluster |>
        dplyr::filter(
            is.na(.data$mz_time_isotope_annotated) |
                .data$mz_time_isotope_annotated != .data$mz_time_annotated
        )

    metabolite_database_simple <- metabolite_database |>
        dplyr::mutate(
            Exact_mass = round(.data$Exact_mass, 8),
            Mono_mass = round(.data$Mono_mass, 4),
            Name = clean_name(.data$Name),
            Most_abundant_isotopologue_mass = round(.data$Most_abundant_isotopologue_mass, 4)
        ) |>
        dplyr::select(
            Exact_mass,
            Mono_mass,
            Name,
            Formula,
            Most_abundant_isotopologue_mass,
            Charge_natural,
            dplyr::any_of(compound_identifier_columns),
            dplyr::any_of(database_identifier_columns),
            Subclass,
            has_NH3,
            has_H2O,
            n_NH3,
            n_H2O
        ) |>
        dplyr::distinct()

    final_results_cluster_all <- final_results_cluster |>
        dplyr::mutate(Exact_mass = round(.data$Exact_mass, 8)) |>
        dplyr::select(-dplyr::any_of("Mono_mass")) |>
        dplyr::inner_join(metabolite_database_simple, by = "Exact_mass") |>
        dplyr::left_join(hmdb_metabolites_concentrations_average_simple, by = "InChIKey") |>
        dplyr::select(
            dplyr::any_of(c(
                "Mono_mass", "identification_type_annotated", "ion_mode",
                "Adduct_annotated", "mz_annotated", "time_annotated",
                "mz_time_annotated", "mean_intensity", "adduct_corr_cluster",
                "mz_isotope", "time_isotope", "mz_time_isotope_annotated",
                "mean_intensity_isotope", "isotopic_correlation",
                "Probability", "cluster_identification"
            )),
            dplyr::any_of(annotation_detail_columns),
            dplyr::any_of(c(
                "Concentration_average", "Concentration_units",
                "has_NH3", "has_H2O", "n_NH3", "n_H2O"
            )),
            dplyr::everything(),
            -dplyr::any_of(exact_mass_columns)
        ) |>
        drop_feature_intensity_columns(feature_intensity_columns)

    write_result_csv(final_results_cluster_all, "clustering_of_adduct_isotope_all_result.csv")

    final_results_cluster_identified <- final_results_cluster |>
        dplyr::filter(
            (.data$Adduct_annotated %in% c("M-H", "M+H") & .data$cluster_identification == 1) |
                .data$Adduct_annotated %in% c("M+Na", "M+K", "M+Cl") |
                .data$Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") |
                .data$Adduct_annotated %in% c("M+2Na-H", "M+2K-H")
        )

    final_results_cluster_identified_2 <- final_results_cluster_identified |>
        dplyr::mutate(Exact_mass = round(.data$Exact_mass, 8)) |>
        dplyr::select(-dplyr::any_of("Mono_mass")) |>
        dplyr::inner_join(metabolite_database_simple, by = "Exact_mass") |>
        dplyr::left_join(hmdb_metabolites_concentrations_average_simple, by = "InChIKey") |>
        dplyr::filter(
            (.data$Adduct_annotated == "M+H-NH3" & .data$has_NH3 == TRUE) |
                (.data$Adduct_annotated %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & .data$n_H2O >= 1) |
                (.data$Adduct_annotated %in% c("M-2H2O+H", "M+H-2H2O") & .data$n_H2O >= 2) |
                (!.data$Adduct_annotated %in% c(
                    "M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H",
                    "M-H-H2O", "M-2H2O+H", "M+H-2H2O"
                ))
        ) |>
        dplyr::filter(
            .data$Adduct_annotated %in% c("M-H", "M+H") |
                (.data$Adduct_annotated %in% c("M+Na", "M+K", "M+Cl") & .data$Subclass == "Carbohydrates and carbohydrate conjugates") |
                (.data$Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & .data$Subclass == "Bile acids, alcohols and derivatives") |
                (.data$Adduct_annotated %in% c("M+H-H2O") & .data$Subclass == "Retinoids") |
                (.data$Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & .data$Subclass == "Pregnane steroids")
        ) |>
        dplyr::rename(
            identification_type = identification_type_annotated,
            Adduct = Adduct_annotated,
            mz = mz_annotated,
            time = time_annotated,
            identified_Name = Name
        ) |>
        dplyr::mutate(
            identification_method = "mz matching; clustering of adducts and isotopes",
            match_category = "multiple"
        ) |>
        dplyr::select(
            dplyr::any_of(c(
                "Mono_mass", "identification_type", "ion_mode", "identified_Name",
                "Formula", "Most_abundant_isotopologue_mass",
                "Charge_natural", "Subclass"
            )),
            dplyr::any_of(c(compound_identifier_columns, database_identifier_columns)),
            dplyr::any_of(c(
                "Adduct", "mz", "time", "cluster_identification",
                "identification_method", "Concentration_average",
                "Concentration_units", "Probability", "match_category"
            )),
            dplyr::everything(),
            -dplyr::any_of(exact_mass_columns)
        ) |>
        dplyr::filter(!is.na(.data$InChIKey)) |>
        dplyr::mutate(
            mz = round(.data$mz, 4),
            time = round(.data$time, 0),
            mz_time = paste0(.data$mz, "_", .data$time)
        ) |>
        dplyr::relocate(mz_time, .after = time)

    write_result_csv(final_results_cluster_identified_2, "clustering_of_adduct_isotope_main_result.csv")

    list(
        mz_only = mz_only_output,
        mz_only_isotope = mz_only_isotope_output,
        clustering_all = final_results_cluster_all,
        clustering_main = final_results_cluster_identified_2
    )
}
