#' @title Get chemical classification for each InChIKey using ClassyFire
#' @description
#'
#' `get_chemical_classification` retrieves chemical classification information for a vector of InChIKey using the ClassyFire API. It processes the classification data and returns a structured data frame with classification levels.
#' 
#' @param unique_inchikey A character vector of InChIKeys to classify.
#' @param query_missing Logical; if `TRUE`, query the ClassyFire API for
#'   InChIKeys not found in the packaged PubChem ClassyFire cache, and for
#'   InChIKeys where cached classification levels are incomplete. Default is
#'   `TRUE`.
#' @param database Deprecated compatibility argument. Cached classification
#'   lookup always uses the packaged PubChem ClassyFire dataset.
#' @param progress_callback Optional function called as
#'   `progress_callback(current, total, detail)` to report classification
#'   progress.
#' @param status_callback Optional function called as
#'   `status_callback(detail)` to report status text.
#' @return A data frame containing chemical classification information with columns: InChIKey, Kingdom, Superclass, Class, Subclass, and Alternative_parent.
#' @export get_chemical_classification

get_chemical_classification = function(unique_inchikey, query_missing = TRUE, database = c("metorigindb", "pubchem"), progress_callback = NULL, status_callback = NULL) {
  library(classyfireR)
  library(dplyr)

  database <- match.arg(database)
  normalize_inchikey <- function(x) {
    toupper(trimws(as.character(x)))
  }

  classification_template = tibble::tibble(
    InChIKey = character(),
    Kingdom = character(),
    Superclass = character(),
    Class = character(),
    Subclass = character(),
    Direct_parent = character(),
    Alternative_parent = character()
  )

  unique_inchikey = normalize_inchikey(unique_inchikey)
  unique_inchikey = unique_inchikey[!is.na(unique_inchikey) & unique_inchikey != ""]
  unique_inchikey = unique(unique_inchikey)

  if (length(unique_inchikey) == 0) {
    return(classification_template)
  }

  report_progress <- function(current, total, detail) {
    if (is.function(progress_callback)) {
      try(progress_callback(current, total, detail), silent = TRUE)
    }
  }

  report_status <- function(detail) {
    if (is.function(status_callback)) {
      try(status_callback(detail), silent = TRUE)
    }
  }

  # create a function to process the classification list
  process_classification_list <- function(classification_list, original_inchikeys) {
    # Initialize an empty data frame to store the results
    ## purrr is needed for the map_df function
    classification_df <- purrr::map2_df(classification_list, original_inchikeys, function(classif, orig_inchikey) {
      # Check if the classification object is NULL
      if (is.null(classif)) {
        # Return original InChIKey with NA values for classification fields
        return(tibble::tibble(
          InChIKey = orig_inchikey,
          Kingdom = NA_character_,
          Superclass = NA_character_,
          Class = NA_character_,
          Subclass = NA_character_,
          Direct_parent = NA_character_,
          Alternative_parent = NA_character_
        ))
      }
      
      # Extract InChIKey from the meta data (but use original as fallback)
      meta_info <- meta(classif)
      inchikey <- meta_info$inchikey
      # Remove "InChIKey=" prefix if present
      inchikey <- sub("InChIKey=", "", inchikey)
      # Use original InChIKey if extraction failed
      if (is.null(inchikey) || is.na(inchikey) || inchikey == "") {
        inchikey <- orig_inchikey
      }
      
      # Extract classification levels
      classif_tibble <- classification(classif)
      
      # Initialize variables with NA
      kingdom <- NA_character_
      superclass <- NA_character_
      class <- NA_character_
      subclass <- NA_character_
      direct_parent <- NA_character_
      alternative_parent <- NA_character_
      # Check and extract each level
      ## kingdom
      if ("kingdom" %in% classif_tibble$Level) {
        kingdom <- classif_tibble$Classification[classif_tibble$Level == "kingdom"]
      }
      ## superclass
      if ("superclass" %in% classif_tibble$Level) {
        superclass <- classif_tibble$Classification[classif_tibble$Level == "superclass"]
      }
      ## class
      if ("class" %in% classif_tibble$Level) {
        class <- classif_tibble$Classification[classif_tibble$Level == "class"]
      }
      ## subclass
      if ("subclass" %in% classif_tibble$Level) {
        subclass <- classif_tibble$Classification[classif_tibble$Level == "subclass"]
      }
      ## direct_parent
      if ("level 7" %in% classif_tibble$Level) {
        direct_parent <- classif_tibble$Classification[classif_tibble$Level == "level 7"]
      } ## if there is no level 7 but there is level 5, then use the level 5 as the direct_parent
      else if (is.na(direct_parent) && "level 5" %in% classif_tibble$Level) {
        direct_parent <- classif_tibble$Classification[classif_tibble$Level == "level 5"]
      } ## if there is no direct_parent, then use the subclass as the direct_parent
      else {
        direct_parent <- subclass
      }
      
      ## alternative_parent (using the alternative_parents function)
      if (length(alternative_parents(classif)) > 0) {
        # extract the name of the alternative_parents
        alternative_parent_name <- alternative_parents(classif)$name
        # merge the alternative_parent_name into a string separated by ";"
        alternative_parent <- paste(alternative_parent_name, collapse = ";")
      }

      # Return a tibble with the extracted information
      tibble::tibble(
        InChIKey = inchikey,
        Kingdom = kingdom,
        Superclass = superclass,
        Class = class,
        Subclass = subclass,
        Direct_parent = direct_parent,
        Alternative_parent = alternative_parent
      )
    })
    
    return(classification_df)
  }

  cached_classification_source <- pubchem_classyfire_loading() |>
    dplyr::transmute(
      InChIKey = .data$InChIKey,
      kingdom = .data$kingdom,
      superclass = .data$superclass,
      class = .data$class,
      subclass = .data$subclass,
      direct_parent = .data$direct_parent,
      alternative_parent = .data$`alternative parents`
    )

  cached_classification_df = cached_classification_source |>
    dplyr::filter(InChIKey %in% unique_inchikey) |>
    dplyr::collect()
  cached_classification_df$InChIKey <- normalize_inchikey(cached_classification_df$InChIKey)

  # keep only classification columns from the cached source
  cached_classification_df_final = cached_classification_df[, c("InChIKey", "kingdom", "superclass", "class", "subclass", "direct_parent", "alternative_parent")]
  # rename the columns
  colnames(cached_classification_df_final) = c("InChIKey", "Kingdom", "Superclass", "Class", "Subclass", "Direct_parent", "Alternative_parent")
  cached_classification_df_final = cached_classification_df_final |>
    dplyr::distinct(InChIKey, .keep_all = TRUE)

  unmatched_inchikey = setdiff(unique_inchikey, cached_classification_df_final$InChIKey)

  normalize_classification_field <- function(x) {
    out <- as.character(x)
    out[trimws(out) == ""] <- NA_character_
    out
  }

  classification_core_columns <- c("Kingdom", "Superclass", "Class", "Subclass")
  cached_classification_df_final <- cached_classification_df_final |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("Kingdom", "Superclass", "Class", "Subclass", "Direct_parent", "Alternative_parent")),
        normalize_classification_field
      )
    )

  incomplete_cache_inchikey <- cached_classification_df_final |>
    dplyr::filter(dplyr::if_any(dplyr::all_of(classification_core_columns), is.na)) |>
    dplyr::pull(InChIKey)

  inchikey_to_query <- unique(c(
    unmatched_inchikey,
    incomplete_cache_inchikey
  ))

  cache_retrieved_count <- length(intersect(unique_inchikey, cached_classification_df_final$InChIKey))
  remaining_to_classify_count <- length(inchikey_to_query)
  progress_total <- max(remaining_to_classify_count + 2L, 2L)

  message("")
  message(
    "PubChem ClassyFire cache retrieved classifications for ",
    cache_retrieved_count,
    " / ",
    length(unique_inchikey),
    " InChIKeys."
  )
  message(
    "Remaining InChIKeys to classify via ClassyFire API: ",
    remaining_to_classify_count,
    "."
  )
  report_status(
    paste0(
      "PubChem ClassyFire cache retrieved classifications for ",
      cache_retrieved_count,
      " / ",
      length(unique_inchikey),
      " InChIKeys. Remaining to classify via ClassyFire API: ",
      remaining_to_classify_count,
      "."
    )
  )
  report_progress(
    1L,
    progress_total,
    paste0(
      "PubChem ClassyFire cache retrieved classifications for ",
      cache_retrieved_count,
      " / ",
      length(unique_inchikey),
      " InChIKeys"
    )
  )

  if (!query_missing || length(inchikey_to_query) == 0) {
    if (!query_missing && remaining_to_classify_count > 0) {
      message("ClassyFire API querying is disabled (`query_missing = FALSE`).")
      report_status("ClassyFire API querying is disabled (`query_missing = FALSE`).")
    }

    missing_keys <- setdiff(unique_inchikey, cached_classification_df_final$InChIKey)
    missing_df = tibble::tibble(
      InChIKey = missing_keys,
      Kingdom = rep(NA_character_, length(missing_keys)),
      Superclass = rep(NA_character_, length(missing_keys)),
      Class = rep(NA_character_, length(missing_keys)),
      Subclass = rep(NA_character_, length(missing_keys)),
      Direct_parent = rep(NA_character_, length(missing_keys)),
      Alternative_parent = rep(NA_character_, length(missing_keys))
    )

    return(
      dplyr::bind_rows(cached_classification_df_final, missing_df) |>
        dplyr::distinct(InChIKey, .keep_all = TRUE)
    )
  }

  n_iter = length(inchikey_to_query)

  message("")
  message("Classifying the chemicals...")
  message("There are ", n_iter, " chemicals to classify.")
  report_status(paste0("Classifying ", n_iter, " chemicals via ClassyFire API..."))
  report_progress(2L, progress_total, paste0("Classifying ", n_iter, " chemicals via ClassyFire API"))

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  pubchem_classification_List <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    current_inchikey = inchikey_to_query[i]
    report_progress(
      2L + i,
      progress_total,
      paste0("Classifying ", i, "/", n_iter, ": ", current_inchikey)
    )
    pubchem_classification_List[[i]] <- tryCatch(
      get_classification(current_inchikey),
      error = function(e) NULL
    )
    utils::setTxtProgressBar(pb, i)
  }

  close(pb)
  report_status("Chemical classification complete.")
  report_progress(progress_total, progress_total, "Chemical classification complete")

  api_classification_df <- process_classification_list(pubchem_classification_List, inchikey_to_query) |>
    dplyr::mutate(
      InChIKey = normalize_inchikey(.data$InChIKey),
      dplyr::across(
        dplyr::all_of(c("Kingdom", "Superclass", "Class", "Subclass", "Direct_parent", "Alternative_parent")),
        normalize_classification_field
      )
    ) |>
    dplyr::distinct(InChIKey, .keep_all = TRUE)

  combined_df <- cached_classification_df_final |>
    dplyr::full_join(api_classification_df, by = "InChIKey", suffix = c("_cache", "_api")) |>
    dplyr::transmute(
      InChIKey = .data$InChIKey,
      Kingdom = dplyr::coalesce(.data$Kingdom_cache, .data$Kingdom_api),
      Superclass = dplyr::coalesce(.data$Superclass_cache, .data$Superclass_api),
      Class = dplyr::coalesce(.data$Class_cache, .data$Class_api),
      Subclass = dplyr::coalesce(.data$Subclass_cache, .data$Subclass_api),
      Direct_parent = dplyr::coalesce(.data$Direct_parent_cache, .data$Direct_parent_api),
      Alternative_parent = dplyr::coalesce(.data$Alternative_parent_cache, .data$Alternative_parent_api)
    ) |>
    dplyr::distinct(InChIKey, .keep_all = TRUE)

  missing_keys <- setdiff(unique_inchikey, combined_df$InChIKey)
  missing_df = tibble::tibble(
    InChIKey = missing_keys,
    Kingdom = rep(NA_character_, length(missing_keys)),
    Superclass = rep(NA_character_, length(missing_keys)),
    Class = rep(NA_character_, length(missing_keys)),
    Subclass = rep(NA_character_, length(missing_keys)),
    Direct_parent = rep(NA_character_, length(missing_keys)),
    Alternative_parent = rep(NA_character_, length(missing_keys))
  )

  dplyr::bind_rows(combined_df, missing_df) |>
    dplyr::distinct(InChIKey, .keep_all = TRUE)
}
