#' @title Filter candidate metabolites by monoisotopic mass
#' @description
#'
#' `mass_filter` filters candidate metabolites from the selected database based
#' on the input m/z, adduct, and monoisotopic mass. The mass-difference
#' estimation window is intentionally wide and estimated from the m/z and ppm
#' threshold. This serves as a preliminary filter before precise m/z matching.
#' 
#' @param mz The m/z value of the features.
#' @param mz_ppm The relative m/z matching threshold for matching features between two data sets. Default is 10 ppm.
#' @param adduct The adduct of the m/z. Default is c("M+H"). If multiple adducts are provided, it would be provided as a vector, like c("M+H", "M+Na").
#' @param database Database used for metabolite lookup. Choose either
#'   `"metorigindb"` or `"pubchem"`. Default is `"metorigindb"`.
#' @return A data frame of candidate metabolites from the selected database.
#' @export mass_filter

mass_filter = function(mz, mz_ppm = 10, adduct = c("M+H"), database = c("metorigindb", "pubchem")){
    library(MetaboCoreUtilsAdduct)
    library(dplyr)
    library(arrow)
    library(tidyr)

    database <- match.arg(database)

    # calculate the theoretical monoisotopic mass given the m/z and adduct
    mass_df = mz2mass_df_safe(x = mz, adduct = adduct)

    ############## estimate the mass difference given the m/z and extract only possible metabolite records from the huge PubChem database
    ## as the mass difference estiamtion window is big, it will retrieve a lot of molecules, but this is a preliminary filter and we will do precise m/z calculation later
    # extract the mz_ppm as specified by the user or the defatul (10 ppm)
    ppm = mz_ppm
    # this is the ppm default for the mass difference estimation
    ppm_default = 10
    # this is the default mass difference estimation window for about 13.33 mz ppm difference given the mass of 1500 Da
    mass_diff_default = 0.02
    # this adjusts the default mass difference estimation window
    new_mass_diff = (ppm/ppm_default)*mass_diff_default
    # print the ppm and new_mass_diff
    message(paste0("Current m/z ppm threshold is: ",ppm))
    message(paste0("Current mass difference estimation window is: ", new_mass_diff, " Da"))

    # calculate the mass_windows using the estimated monoisotopic mass from mass_df and the mass difference estimation window
    mass_window_low = mass_df$adduct_mass - new_mass_diff
    mass_window_high = mass_df$adduct_mass + new_mass_diff
    mass_window = data.frame(low = mass_window_low, high = mass_window_high)

    # in case where the monoisotopic mass is not neutral, we can minus the 16 masses of proton (1.007276 Da) from the mass_window (the maximum absolute number of natural charges is 16)
    mass_window_low_neutral = mass_window_low - 16*1.007276
    mass_window_high_neutral = mass_window_high - 16*1.007276
    mass_window_neutral = data.frame(low = mass_window_low_neutral, high = mass_window_high_neutral)

    # combine the mass_window and mass_window_neutral
    mass_window_combined = rbind(mass_window, mass_window_neutral)

    # build one OR-ed expression string for each pair of mass window (basically each row)
    pred <- purrr::map2(
      mass_window_combined$low, mass_window_combined$high,
      ~ expr(between(Mono_mass, !!.x, !!.y))
    ) |>
      purrr::reduce(~ expr((!!.x) | (!!.y)))

    # pull the molecules out of the selected dataset based on the estimated mass window
    mass_candidate = metabolite_database_loading(database = database) |>
      filter(!!pred) |>            # Arrow evaluates it, no RAM blow-up
      dplyr::collect()

    # cross join Adduct with mass_candidate to form a adduct list for each unique KEGGID 
    ## if there is only one adduct, then only one adduct will be added to the mass_candidate
    ## if there are multiple adducts, then each adduct will be added to the mass_candidate so that each metabolite will have multiple adducts
    df_adduct <- data.frame(Adduct = adduct, stringsAsFactors = FALSE)

    # do the Cartesian join
    mass_candidate_2 <- merge(
        df_adduct,
        mass_candidate,
        by = NULL
    )

    # convert mass_candidate_2 to a tibble
    mass_candidate_2 = tibble(mass_candidate_2)

    # return the mass candidate with the adduct list
    return(mass_candidate_2)
}
