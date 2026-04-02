#' @title Run Lightweight m/z Clustering with PubChem
#' @description
#' `mz_match_clustering_pubchem()` is a PubChem-focused variant of
#' `mz_match_clustering()` designed to reduce memory and startup time.
#'
#' It applies three lightweight steps before clustering:
#' 1) Arrow-side neutral-mass range filtering before `collect()`.
#' 2) Adduct-wise streaming pre-screening (no full PubChem x adduct expansion).
#' 3) Minimal-column loading for clustering-required fields.
#'
#' After pre-screening, it calls `mz_match_clustering()` using the reduced
#' candidate table as `metabolite_database`.
#'
#' @inheritParams mz_match_clustering
#' @return Same return value as `mz_match_clustering()`.
#' @export
mz_match_clustering_pubchem <- function(
    met_raw_wide,
    metabolite_database = NULL,
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

    default_adducts <- function(current_ion_mode) {
        if (identical(current_ion_mode, "negative")) {
            return(c("M-H", "M+Cl", "M-H-H2O", "2M-H"))
        }

        c(
            "M+H", "M+Na", "M+2Na-H", "M+H-H2O", "M+H-NH3",
            "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O"
        )
    }

    adduct_allowed_by_formula <- function(data, adduct) {
        if (identical(adduct, "M+H-NH3")) {
            return(data$has_NH3 %in% TRUE)
        }
        if (adduct %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O")) {
            return(data$n_H2O >= 1)
        }
        if (adduct %in% c("M-2H2O+H", "M+H-2H2O")) {
            return(data$n_H2O >= 2)
        }

        rep(TRUE, nrow(data))
    }

    if (is.null(All_Adduct)) {
        All_Adduct <- default_adducts(ion_mode)
    }

    if (!is.null(metabolite_database)) {
        if (isTRUE(show_progress)) {
            message("Using provided PubChem metabolite_database (skipping lightweight prefilter).")
        }

        return(mz_match_clustering(
            met_raw_wide = met_raw_wide,
            metabolite_database = metabolite_database,
            database = "pubchem",
            hmdb_concentration = hmdb_concentration,
            mz_threshold = mz_threshold,
            biospecimen = biospecimen,
            All_Adduct = All_Adduct,
            ion_mode = ion_mode,
            adduct_correlation_r_threshold = adduct_correlation_r_threshold,
            adduct_correlation_time_threshold = adduct_correlation_time_threshold,
            isotopic_correlation_r_threshold = isotopic_correlation_r_threshold,
            isotopic_correlation_time_threshold = isotopic_correlation_time_threshold,
            imputation_method = imputation_method,
            output_dir = output_dir,
            write_output = write_output,
            show_progress = show_progress,
            progress_callback = progress_callback
        ))
    }

    met_raw_wide <- dplyr::as_tibble(met_raw_wide)
    if (ncol(met_raw_wide) < 2) {
        stop("met_raw_wide must contain at least mz and time columns.", call. = FALSE)
    }
    colnames(met_raw_wide)[1:2] <- c("mz", "time")

    mz_values <- suppressWarnings(as.numeric(met_raw_wide$mz))
    mz_values <- mz_values[is.finite(mz_values)]
    if (length(mz_values) == 0) {
        stop("No finite mz values found in met_raw_wide.", call. = FALSE)
    }

    adduct_definition <- adduct_definition_loading() |>
        dplyr::filter(.data$name %in% All_Adduct) |>
        dplyr::select(name, mass_multi, mass_add)

    if (nrow(adduct_definition) == 0) {
        stop("None of the selected adducts are available in adduct_definition.", call. = FALSE)
    }

    ppm_margin <- max(mz_threshold, 0) / 1000000
    mz_lower <- min(mz_values) * (1 - ppm_margin)
    mz_upper <- max(mz_values) * (1 + ppm_margin)

    mass_at_mz_lower <- (mz_lower - adduct_definition$mass_add) / adduct_definition$mass_multi
    mass_at_mz_upper <- (mz_upper - adduct_definition$mass_add) / adduct_definition$mass_multi
    lower_candidates <- pmin(mass_at_mz_lower, mass_at_mz_upper)
    upper_candidates <- pmax(mass_at_mz_lower, mass_at_mz_upper)
    neutral_mass_lower <- min(pmin(lower_candidates, upper_candidates), na.rm = TRUE)
    neutral_mass_upper <- max(pmax(lower_candidates, upper_candidates), na.rm = TRUE)

    # include isotopologue range conservatively
    neutral_mass_lower <- neutral_mass_lower - 1.01
    neutral_mass_upper <- neutral_mass_upper + 1.01

    cache_key <- paste0(
        "pubchem_clustering_prefilter::",
        round(neutral_mass_lower, 6), "::",
        round(neutral_mass_upper, 6), "::",
        paste(sort(All_Adduct), collapse = "|")
    )
    pubchem_prefiltered <- massmatcher_cache_get(cache_key)

    if (is.null(pubchem_prefiltered)) {
        if (isTRUE(show_progress)) {
            message("Prefiltering PubChem by neutral mass range before collect...")
        }

        location <- locate_database_parquet("pubchem")
        dataset_key <- paste0(
            "metabolite_raw_dataset::pubchem::",
            normalizePath(location$path, winslash = "/", mustWork = FALSE)
        )
        dataset <- massmatcher_cache_get(dataset_key)
        if (is.null(dataset)) {
            message(location$message)
            dataset <- arrow::open_dataset(location$path)
            massmatcher_cache_set(dataset_key, dataset)
        }

        pubchem_prefiltered <- dataset |>
            dplyr::transmute(
                Name = .data$Name,
                Formula = .data$Formula,
                Mono_mass = .data$MonoMass,
                Most_abundant_isotopologue_mass = .data$most_abundant_isotopologue_mass,
                InChIKey = .data$InChIKey,
                CID = .data$CID,
                HMDB_ID = .data$HMDB_ID,
                KEGG_ID = .data$KEGG_ID
            ) |>
            dplyr::filter(
                (.data$Mono_mass >= neutral_mass_lower & .data$Mono_mass <= neutral_mass_upper) |
                    (
                        .data$Most_abundant_isotopologue_mass >= neutral_mass_lower &
                            .data$Most_abundant_isotopologue_mass <= neutral_mass_upper
                    )
            ) |>
            dplyr::collect() |>
            dplyr::mutate(
                Subclass = NA_character_,
                has_NH3 = FALSE,
                n_NH3 = 0L,
                has_H2O = FALSE,
                n_H2O = 0L,
                Mass_diff = NA_real_
            ) |>
            dplyr::filter(!is.na(.data$Mono_mass), is.finite(.data$Mono_mass)) |>
            dplyr::distinct(.data$InChIKey, .keep_all = TRUE)

        massmatcher_cache_set(cache_key, pubchem_prefiltered)
    } else if (isTRUE(show_progress)) {
        message("Using cached PubChem prefilter for this mass/adduct window...")
    }

    if (nrow(pubchem_prefiltered) > 0) {
        if (isTRUE(show_progress)) {
            message("Running adduct-wise streaming pre-screen to avoid full PubChem x adduct expansion...")
        }

        pubchem_prefiltered <- dplyr::mutate(pubchem_prefiltered, .row_id = dplyr::row_number())
        matched_row_ids <- integer(0)

        for (adduct_name in All_Adduct) {
            keep_idx <- adduct_allowed_by_formula(pubchem_prefiltered, adduct_name)
            if (!any(keep_idx)) {
                next
            }

            db_sub <- pubchem_prefiltered[keep_idx, , drop = FALSE]

            mono_mz <- mass2mz_df_safe(mass = db_sub$Mono_mass, adduct = adduct_name)$mz
            mono_keep <- which(is.finite(mono_mz))
            if (length(mono_keep) > 0) {
                mono_overlap <- find.Overlapping.mzs(
                    data.frame(mz = mz_values, stringsAsFactors = FALSE),
                    data.frame(mz = mono_mz[mono_keep], stringsAsFactors = FALSE),
                    mz.thresh = mz_threshold
                )
                if (nrow(mono_overlap) > 0) {
                    matched_row_ids <- c(matched_row_ids, db_sub$.row_id[mono_keep][mono_overlap$index.B])
                }
            }

            if (any(is.finite(db_sub$Most_abundant_isotopologue_mass))) {
                iso_mz <- mass2mz_df_safe(
                    mass = db_sub$Most_abundant_isotopologue_mass,
                    adduct = adduct_name
                )$mz
                iso_keep <- which(is.finite(iso_mz))
                if (length(iso_keep) > 0) {
                    iso_overlap <- find.Overlapping.mzs(
                        data.frame(mz = mz_values, stringsAsFactors = FALSE),
                        data.frame(mz = iso_mz[iso_keep], stringsAsFactors = FALSE),
                        mz.thresh = mz_threshold
                    )
                    if (nrow(iso_overlap) > 0) {
                        matched_row_ids <- c(matched_row_ids, db_sub$.row_id[iso_keep][iso_overlap$index.B])
                    }
                }
            }
        }

        matched_row_ids <- unique(matched_row_ids)
        pubchem_prefiltered <- pubchem_prefiltered |>
            dplyr::filter(.data$.row_id %in% matched_row_ids) |>
            dplyr::select(-.data$.row_id)
    }

    if (isTRUE(show_progress)) {
        message("Lightweight PubChem candidates prepared: ", nrow(pubchem_prefiltered), " rows.")
    }

    mz_match_clustering(
        met_raw_wide = met_raw_wide,
        metabolite_database = pubchem_prefiltered,
        database = "pubchem",
        hmdb_concentration = hmdb_concentration,
        mz_threshold = mz_threshold,
        biospecimen = biospecimen,
        All_Adduct = All_Adduct,
        ion_mode = ion_mode,
        adduct_correlation_r_threshold = adduct_correlation_r_threshold,
        adduct_correlation_time_threshold = adduct_correlation_time_threshold,
        isotopic_correlation_r_threshold = isotopic_correlation_r_threshold,
        isotopic_correlation_time_threshold = isotopic_correlation_time_threshold,
        imputation_method = imputation_method,
        output_dir = output_dir,
        write_output = write_output,
        show_progress = show_progress,
        progress_callback = progress_callback
    )
}
