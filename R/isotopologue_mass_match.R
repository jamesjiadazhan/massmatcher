#' @title Match m/z values to candidate metabolites by isotopologue mass
#' @description
#'
#' `isotopologue_mass_match` first applies the same broad exact-mass database
#' prefiltering strategy used by `mz_match_clustering()` when
#' `metabolite_database` is not supplied, then performs exact isotopologue m/z
#' matching within the specified ppm threshold.
#' 
#' @param unknown_feature The mz and time of the unknown features. The first column name has to be  mz and the second column name has to be time.
#' @param mz_ppm The m/z ppm threshold for matching the m/z values to the candidate metabolites. Default is 10 ppm.
#' @param adduct The adduct of the m/z. Default is c("M+H"). If multiple adducts are provided, it would be provided as a vector, like c("M+H", "M+Na").
#' @param database Database used for metabolite lookup. Choose either
#'   `"metorigindb"` or `"pubchem"`. Default is `"metorigindb"`.
#' @param metabolite_database Optional prefiltered metabolite database. When
#'   provided, `isotopologue_mass_match()` uses it instead of querying the full
#'   selected database.
#' @return A data frame of mass matched metabolites from the selected database.
#' @export isotopologue_mass_match

isotopologue_mass_match = function(unknown_feature, mz_ppm = 10, adduct = c("M+H"), database = c("metorigindb", "pubchem"), metabolite_database = NULL){
    library(MetaboCoreUtilsAdduct)
    library(dplyr)

    database <- match.arg(database)

    if (is.null(metabolite_database)) {
        metabolite_database <- prefilter_metabolite_database_by_mz_range(
            mz_values = unknown_feature$mz,
            adducts = adduct,
            ppm_threshold = mz_ppm,
            database = database
        )
    }
    
    mass_candidate <- tidyr::crossing(
        dplyr::as_tibble(metabolite_database),
        Adduct = adduct
    )

    # extract the mz from the unknown_feature
    mz = unknown_feature$mz

    ########### Calculate the theoretical m/z values for the mass candidate given the adduct
    mass_candidate_theoretical_mz = mass2mz_df_safe(
        mass = mass_candidate$Exact_mass_most_abundant_isotopologue,
        adduct = mass_candidate$Adduct
    )
    mass_candidate_combined = mass_candidate
    mass_candidate_combined$theoretical_mz = mass_candidate_theoretical_mz$mz

    ########### mz matching for between the given m/z and the mass candidate from PubChem database
    # extract the mz from the unknown_feature
    mz = unknown_feature$mz
    # select the mz column
    mass_candidate_mz = mass_candidate_combined[,"theoretical_mz"]
    # find the overlapping mz within the specified m/z ppm threshold
    masteroverlap.mass_candidate = find.Overlapping.mzs(mz, mass_candidate_mz, mz.thresh = mz_ppm)
    # select the matched unknown_feature and mass candidate from the mass candidate data frame
    mz_matched = dplyr::slice(unknown_feature, masteroverlap.mass_candidate$index.A)
    # select the matched mass candidate from the mass candidate data frame
    mass_candidate_matched = dplyr::slice(mass_candidate_combined, masteroverlap.mass_candidate$index.B)
    # combine the mz_matched and mass_candidate_matched into a data frame
    mass_candidate_matched_final = cbind(mz_matched, mass_candidate_matched)
    # convert the mass_candidate_matched_final to a tibble
    mass_candidate_matched_final = tibble::as_tibble(mass_candidate_matched_final)
    # relocate the theoretical_mz column right after the time column
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::relocate(theoretical_mz, .after = time)
    # calculate the absolute mz ppm difference between the mz and theoretical_mz
    mass_candidate_matched_final$mz_ppm = abs((mass_candidate_matched_final$mz - mass_candidate_matched_final$theoretical_mz) / mass_candidate_matched_final$theoretical_mz * 1000000)
    # round the mz_ppm to 2 decimal places
    mass_candidate_matched_final$mz_ppm = round(mass_candidate_matched_final$mz_ppm, 2)
    # relocate the mz_ppm column to the third column
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::relocate(mz_ppm, .after = theoretical_mz)
    # keep Adduct close to the ppm error for easier scanning in UI/export
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::relocate(dplyr::any_of("Adduct"), .after = mz_ppm)
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::select(-dplyr::any_of(c("Exact_mass", "Exact_mass_most_abundant_isotopologue", "Mass_diff")))
    # arrange the result by time and then mz
    mass_candidate_matched_final = mass_candidate_matched_final %>%
        dplyr::arrange(time) %>%
        dplyr::arrange(mz)
    # return the mass_matched data frame
    return(mass_candidate_matched_final)
}
