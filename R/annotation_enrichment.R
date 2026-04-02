#' @title Add chemical classification and concentration to annotation results
#' @description
#' `annotate_match_results` enriches result tables (for example outputs from
#' `mass_match`, `isotopologue_mass_match`, or `mz_match_clustering`) by adding
#' ClassyFire-style classification columns and HMDB concentration columns in one
#' step.
#'
#' @param result_table A data frame/tibble containing at least one of
#'   `InChIKey` or `HMDB_ID`.
#' @param database Database used for lookup. Choose either `"metorigindb"` or
#'   `"pubchem"`. Default is `"metorigindb"`.
#' @param hmdb_concentration Optional HMDB concentration table. If `NULL`,
#'   `hmdb_concentration_loading()` is used.
#' @param biospecimen Biospecimen used to prioritize HMDB concentrations when
#'   the concentration table contains a `Biospecimen` column. Default is
#'   `"Blood"`.
#' @param query_missing_classification Logical; if `TRUE`, query ClassyFire API
#'   for unmatched or incompletely-classified InChIKeys. Default is `TRUE`.
#' @param progress_callback Optional function called as
#'   `progress_callback(current, total, detail)` to report classification
#'   progress to an external UI.
#' @param status_callback Optional function called as `status_callback(detail)`
#'   to report enrichment status text.
#' @return An enriched tibble with classification and concentration columns
#'   added when possible.
#' @export
annotate_match_results <- function(
    result_table,
    database = c("metorigindb", "pubchem"),
    hmdb_concentration = NULL,
    biospecimen = "Blood",
    query_missing_classification = TRUE,
    progress_callback = NULL,
    status_callback = NULL
) {
    library(dplyr)

    database <- match.arg(database)
    result <- tibble::as_tibble(result_table)
    classification_columns <- c(
        "Classification_Kingdom",
        "Classification_Superclass",
        "Classification_Class",
        "Classification_Subclass",
        "Classification_Direct_parent",
        "Classification_Alternative_parent"
    )

    if (nrow(result) == 0) {
        for (classification_col in classification_columns) {
            if (!classification_col %in% colnames(result)) {
                result[[classification_col]] <- character(0)
            }
        }
        if (!"Concentration_average" %in% colnames(result)) {
            result$Concentration_average <- numeric(0)
        }
        if (!"Concentration_units" %in% colnames(result)) {
            result$Concentration_units <- character(0)
        }
        return(result)
    }

    if (!"InChIKey" %in% colnames(result) && "INCHIKEY_ID" %in% colnames(result)) {
        result$InChIKey <- result$INCHIKEY_ID
    }
    if (!"HMDB_ID" %in% colnames(result) && "HMDBID" %in% colnames(result)) {
        result$HMDB_ID <- result$HMDBID
    }

    has_inchikey <- "InChIKey" %in% colnames(result)
    has_hmdb <- "HMDB_ID" %in% colnames(result)

    if (!has_inchikey && !has_hmdb) {
        return(result)
    }

    # normalize InChIKey for exact cache lookup
    if (has_inchikey) {
        result <- result |>
            dplyr::mutate(
                .InChIKey_norm = toupper(trimws(as.character(.data$InChIKey)))
            )
    } else {
        result$.InChIKey_norm <- NA_character_
    }

    # Classification enrichment
    unique_inchikey <- unique(result$.InChIKey_norm[!is.na(result$.InChIKey_norm) & result$.InChIKey_norm != ""])
    if (length(unique_inchikey) > 0) {
        class_df <- get_chemical_classification(
            unique_inchikey = unique_inchikey,
            query_missing = query_missing_classification,
            database = database,
            progress_callback = progress_callback,
            status_callback = status_callback
        ) |>
            dplyr::mutate(.InChIKey_norm = toupper(trimws(as.character(.data$InChIKey)))) |>
            dplyr::select(
                .InChIKey_norm,
                Classification_Kingdom = Kingdom,
                Classification_Superclass = Superclass,
                Classification_Class = Class,
                Classification_Subclass = Subclass,
                Classification_Direct_parent = Direct_parent,
                Classification_Alternative_parent = Alternative_parent
            ) |>
            dplyr::distinct(.data$.InChIKey_norm, .keep_all = TRUE)

        result <- result |>
            dplyr::left_join(class_df, by = ".InChIKey_norm")
    }

    for (classification_col in classification_columns) {
        if (!classification_col %in% colnames(result)) {
            result[[classification_col]] <- NA_character_
        }
    }

    # Concentration enrichment
    concentration_df <- hmdb_concentration
    if (is.null(concentration_df)) {
        concentration_df <- hmdb_concentration_loading()
    }
    concentration_df <- dplyr::as_tibble(concentration_df)

    hmdb_source_col <- if ("HMDBID" %in% colnames(concentration_df)) {
        "HMDBID"
    } else if ("HMDB_ID" %in% colnames(concentration_df)) {
        "HMDB_ID"
    } else {
        NA_character_
    }

    if (!is.na(hmdb_source_col) &&
        all(c(hmdb_source_col, "Concentration_average", "Concentration_units") %in% colnames(concentration_df))) {

        if ("Biospecimen" %in% colnames(concentration_df)) {
            concentration_df <- concentration_df |>
                dplyr::filter(.data$Biospecimen == biospecimen)
        }

        concentration_lookup <- concentration_df |>
            dplyr::select(
                HMDB_ID_lookup = dplyr::all_of(hmdb_source_col),
                Concentration_average_lookup = Concentration_average,
                Concentration_units_lookup = Concentration_units
            ) |>
            dplyr::filter(!is.na(.data$HMDB_ID_lookup)) |>
            dplyr::group_by(.data$HMDB_ID_lookup) |>
            dplyr::slice_max(.data$Concentration_average_lookup, n = 1, with_ties = FALSE) |>
            dplyr::ungroup()

        result$HMDB_ID_lookup <- if ("HMDB_ID" %in% colnames(result)) {
            as.character(result$HMDB_ID)
        } else {
            NA_character_
        }

        # Fill missing HMDB_ID from database mapping via exact InChIKey if needed.
        if (any(is.na(result$HMDB_ID_lookup) | result$HMDB_ID_lookup == "")) {
            db_hmdb_map <- metabolite_database_loading(database = database) |>
                dplyr::select(
                    InChIKey,
                    HMDB_ID_map = HMDB_ID
                ) |>
                dplyr::collect() |>
                dplyr::mutate(
                    .InChIKey_norm = toupper(trimws(as.character(.data$InChIKey))),
                    HMDB_ID_map = as.character(.data$HMDB_ID_map)
                ) |>
                dplyr::filter(!is.na(.data$HMDB_ID_map), .data$HMDB_ID_map != "", !is.na(.data$.InChIKey_norm), .data$.InChIKey_norm != "") |>
                dplyr::select(.InChIKey_norm, HMDB_ID_map) |>
                dplyr::distinct(.data$.InChIKey_norm, .keep_all = TRUE)

            result <- result |>
                dplyr::left_join(db_hmdb_map, by = ".InChIKey_norm") |>
                dplyr::mutate(HMDB_ID_lookup = dplyr::coalesce(.data$HMDB_ID_lookup, .data$HMDB_ID_map))
        }

        result <- result |>
            dplyr::left_join(concentration_lookup, by = "HMDB_ID_lookup")

        if ("Concentration_average" %in% colnames(result)) {
            result$Concentration_average <- dplyr::coalesce(
                suppressWarnings(as.numeric(result$Concentration_average)),
                suppressWarnings(as.numeric(result$Concentration_average_lookup))
            )
        } else {
            result$Concentration_average <- suppressWarnings(as.numeric(result$Concentration_average_lookup))
        }

        if ("Concentration_units" %in% colnames(result)) {
            result$Concentration_units <- dplyr::coalesce(
                as.character(result$Concentration_units),
                as.character(result$Concentration_units_lookup)
            )
        } else {
            result$Concentration_units <- as.character(result$Concentration_units_lookup)
        }
    }

    result |>
        dplyr::select(
            -dplyr::any_of(c(
                ".InChIKey_norm",
                "HMDB_ID_map",
                "HMDB_ID_lookup",
                "Concentration_average_lookup",
                "Concentration_units_lookup"
            ))
        )
}

#' @title Add chemical classification and concentration to clustering outputs
#' @description
#' `annotate_mz_match_clustering_results` enriches existing
#' `mz_match_clustering()` output tables (`mz_only`, `mz_only_isotope`,
#' `clustering_all`, and `clustering_main`) with classification and
#' concentration columns via `annotate_match_results`, without rerunning
#' clustering.
#'
#' @param clustering_output A named list returned by `mz_match_clustering()`.
#' @inheritParams annotate_match_results
#' @return A named list with enriched clustering result tables.
#' @export
annotate_mz_match_clustering_results <- function(
    clustering_output,
    database = c("metorigindb", "pubchem"),
    hmdb_concentration = NULL,
    biospecimen = "Blood",
    query_missing_classification = TRUE
) {
    database <- match.arg(database)

    if (!is.list(clustering_output)) {
        stop("`clustering_output` must be a named list from `mz_match_clustering()`.")
    }

    out <- clustering_output
    expected_tables <- c("mz_only", "mz_only_isotope", "clustering_all", "clustering_main")
    missing_tables <- setdiff(expected_tables, names(out))
    if (length(missing_tables) > 0) {
        warning(
            "Missing expected clustering output tables: ",
            paste(missing_tables, collapse = ", "),
            ". Available tables will still be enriched."
        )
    }

    for (nm in intersect(expected_tables, names(out))) {
        out[[nm]] <- annotate_match_results(
            result_table = out[[nm]],
            database = database,
            hmdb_concentration = hmdb_concentration,
            biospecimen = biospecimen,
            query_missing_classification = query_missing_classification
        )
    }

    out
}
