#' Prefilter metabolite database by a broad m/z-derived neutral-mass range
#'
#' Internal/search helper that mirrors the coarse database prefiltering strategy
#' used by `mz_match_clustering()`: determine one broad neutral-mass interval
#' from all submitted m/z values and adducts, filter the Arrow dataset once, and
#' collect only that reduced subset into R.
#'
#' @param mz_values Numeric vector of m/z values.
#' @param adducts Character vector of adducts.
#' @param ppm_threshold Numeric ppm tolerance.
#' @param database Database used for metabolite lookup. Choose either
#'   `"metorigindb"` or `"pubchem"`.
#' @return A tibble/data frame containing a broad candidate subset.
#' @export
prefilter_metabolite_database_by_mz_range <- function(
    mz_values,
    adducts,
    ppm_threshold = 10,
    database = c("metorigindb", "pubchem")
) {
    database <- match.arg(database)

    mz_values <- suppressWarnings(as.numeric(mz_values))
    mz_values <- mz_values[is.finite(mz_values)]
    if (length(mz_values) == 0) {
        stop("`mz_values` must contain at least one finite numeric m/z value.", call. = FALSE)
    }

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

    metabolite_database_loading(database = database) |>
        dplyr::filter(
            (.data$Exact_mass >= lower & .data$Exact_mass <= upper) |
                (
                    .data$Exact_mass_most_abundant_isotopologue >= lower &
                        .data$Exact_mass_most_abundant_isotopologue <= upper
                )
        ) |>
        dplyr::collect()
}
