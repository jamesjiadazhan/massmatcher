#' Merge overlapping numeric intervals
#'
#' Internal helper used to keep Arrow-side mass-window predicates small enough
#' to stay compatible with more restrictive Arrow backends used in deployment.
merge_mass_windows <- function(low, high) {
  windows <- data.frame(
    low = suppressWarnings(as.numeric(low)),
    high = suppressWarnings(as.numeric(high)),
    stringsAsFactors = FALSE
  )
  windows <- windows[is.finite(windows$low) & is.finite(windows$high), , drop = FALSE]
  if (nrow(windows) == 0) {
    return(windows)
  }

  lower <- pmin(windows$low, windows$high)
  upper <- pmax(windows$low, windows$high)
  windows <- windows[order(lower, upper), , drop = FALSE]
  windows$low <- lower[order(lower, upper)]
  windows$high <- upper[order(lower, upper)]

  merged <- vector("list", nrow(windows))
  current_low <- windows$low[[1]]
  current_high <- windows$high[[1]]
  out_idx <- 1L

  if (nrow(windows) > 1) {
    for (i in 2:nrow(windows)) {
      next_low <- windows$low[[i]]
      next_high <- windows$high[[i]]

      if (next_low <= current_high) {
        current_high <- max(current_high, next_high)
      } else {
        merged[[out_idx]] <- data.frame(low = current_low, high = current_high)
        out_idx <- out_idx + 1L
        current_low <- next_low
        current_high <- next_high
      }
    }
  }

  merged[[out_idx]] <- data.frame(low = current_low, high = current_high)
  merged <- merged[seq_len(out_idx)]
  do.call(rbind, merged)
}

#' Collect dataset rows that fall within merged mass windows
#'
#' Internal helper to avoid one very large Arrow predicate when users submit
#' multiple m/z values in the Shiny search module.
collect_candidates_by_mass_windows <- function(dataset, mass_column, low, high, max_windows_per_query = 4L) {
  windows <- merge_mass_windows(low = low, high = high)
  if (nrow(windows) == 0) {
    return(dplyr::collect(dplyr::slice_head(dataset, n = 0)))
  }

  batch_ids <- ceiling(seq_len(nrow(windows)) / max_windows_per_query)
  collected <- vector("list", max(batch_ids))

  for (batch_id in seq_len(max(batch_ids))) {
    batch_windows <- windows[batch_ids == batch_id, , drop = FALSE]
    pred <- purrr::map2(
      batch_windows$low,
      batch_windows$high,
      ~ rlang::expr(dplyr::between(!!rlang::sym(mass_column), !!.x, !!.y))
    ) |>
      purrr::reduce(~ rlang::expr((!!.x) | (!!.y)))

    collected[[batch_id]] <- dataset |>
      dplyr::filter(!!pred) |>
      dplyr::collect()
  }

  collected <- collected[vapply(collected, nrow, integer(1)) > 0]
  if (length(collected) == 0) {
    return(dplyr::collect(dplyr::slice_head(dataset, n = 0)))
  }

  dplyr::bind_rows(collected) |>
    dplyr::distinct()
}
