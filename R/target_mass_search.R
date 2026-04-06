#' @title Target mass search for the given mass and adduct
#' @description
#'
#' `target_mass_search` is a function that searches for the given mass and adduct in the feature table.
#' 
#' @param mass The mass of the targeted chemical.
#' @param feature_table The feature table of the untargeted features, assuming the m/z value column is named "mz".
#' @param adduct The hypothesized adduct of the chemical.
#' @param mz_ppm The m/z ppm threshold for matching the m/z values to the candidate metabolites. Default is 10 ppm.
#' @return A data frame of mass matched (within the specified m/z ppm threshold) features from the feature table.
#' @export
target_mass_search <- function(mass, feature_table, adduct = c("M+H"), mz_ppm = 10) {
    library(MetaboCoreUtilsAdduct)
    library(dplyr)
    library(tidyr)

    mass <- suppressWarnings(as.numeric(mass))
    mass <- mass[is.finite(mass)]
    if (length(mass) == 0) {
        stop("`mass` must include at least one numeric value.", call. = FALSE)
    }

    feature_table <- as.data.frame(feature_table, stringsAsFactors = FALSE)
    if (!("mz" %in% names(feature_table))) {
        stop("`feature_table` must include an `mz` column.", call. = FALSE)
    }
    feature_table$mz <- suppressWarnings(as.numeric(feature_table$mz))
    if ("time" %in% names(feature_table)) {
        feature_table$time <- suppressWarnings(as.numeric(feature_table$time))
    }
    feature_table <- feature_table[is.finite(feature_table$mz), , drop = FALSE]
    if (nrow(feature_table) == 0) {
        stop("`feature_table$mz` has no valid numeric values.", call. = FALSE)
    }

    # crossing the mass and adduct
    mass_candidate = tidyr::crossing(mass, adduct)

    # calculate the theoretical m/z values 
    mass_candidate_mz = mass2mz_df_safe(mass = mass_candidate$mass, adduct = mass_candidate$adduct)

    # add the theoretical m/z values to the mass_candidate
    mass_candidate$theoretical_mz = mass_candidate_mz$mz

    ########### mz matching for between the given m/z and the mass candidate from the feature table
    # select the mz column from the feature table
    mz = feature_table$mz
    # select the theoretical m/z column
    mass_candidate_mz = mass_candidate$theoretical_mz
    # find the overlapping mz within the specified m/z ppm threshold
    masteroverlap.mass_candidate = find.Overlapping.mzs(mz, mass_candidate_mz, mz.thresh = mz_ppm)
    # select the matched mz in the feature table
    mz_matched = dplyr::slice(feature_table, masteroverlap.mass_candidate$index.A)
    # select the matched mass candidate
    mass_candidate_matched = dplyr::slice(mass_candidate, masteroverlap.mass_candidate$index.B)
    # combine the mz_matched and mass_candidate_matched into a data frame
    mass_candidate_matched_final = cbind(mass_candidate_matched, mz_matched)
    # convert the mass_candidate_matched_final to a tibble
    mass_candidate_matched_final = tibble::as_tibble(mass_candidate_matched_final)
    # calculate the absolute mz ppm difference between the mz and theoretical_mz
    mass_candidate_matched_final$mz_ppm = abs((mass_candidate_matched_final$mz - mass_candidate_matched_final$theoretical_mz) / mass_candidate_matched_final$theoretical_mz * 1000000)
    # round the mz_ppm to 2 decimal places
    mass_candidate_matched_final$mz_ppm = round(mass_candidate_matched_final$mz_ppm, 2)
    # relocate the mz_ppm column to the third column
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::relocate(mz_ppm, .after = theoretical_mz)
    # return the mass_matched data frame
    return(mass_candidate_matched_final)
}
