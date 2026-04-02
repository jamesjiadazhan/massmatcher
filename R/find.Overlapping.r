#' find.Overlapping.mzs
#' 
#' This function finds overlapping m/z and/or retention time values between two data sets.
#' Originally written by Dr. Karan Uppal, Emory University.
#' Enhanced by James Zhan using data.table.
#' @import data.table
#' @param dataA A feature table. The first column is considered as mz. If the second column
#'   is included (optional), it is considered as retention time.
#' @param dataB A feature table. The first column is considered as mz. If the second column
#'   is included (optional), it is considered as retention time.
#' @param mz.thresh The m/z threshold for matching features between two data sets.
#'   Default is 5 ppm.
#' @param time.thresh The retention time threshold for matching features between two data sets.
#'   Default is NA (retention time is not used for matching). If this is used,
#'   the recommended input is 30 (seconds).
#' @return Matrix of overlapping features with columns:
#'   index.A, mz.data.A, (optional) time.data.A, index.B, mz.data.B,
#'   (optional) time.data.B, (optional) time.difference.
#' @export find.Overlapping.mzs

find.Overlapping.mzs <- function(dataA, dataB, mz.thresh = 5, time.thresh = NA) {
    to_feature_dt <- function(x, label) {
        if (data.table::is.data.table(x)) {
            dt <- data.table::copy(x)
        } else if (is.data.frame(x)) {
            dt <- data.table::as.data.table(x)
        } else if (is.matrix(x)) {
            dt <- data.table::as.data.table(as.data.frame(x, stringsAsFactors = FALSE))
        } else if (is.atomic(x) || is.list(x)) {
            dt <- data.table::data.table(mz = as.vector(x))
        } else {
            stop(label, " must be a vector, matrix, data.frame, or data.table.", call. = FALSE)
        }

        if (ncol(dt) < 1) {
            stop(label, " must contain at least one column (m/z).", call. = FALSE)
        }

        dt[, .original_index__ := .I]
        dt
    }

    choose_numeric_column <- function(dt, label, kind = c("m/z", "time"), excluded = character(0)) {
        kind <- match.arg(kind)
        available_cols <- setdiff(names(dt), c(".original_index__", excluded))
        if (length(available_cols) == 0) {
            stop(label, " has no usable columns for ", kind, ".", call. = FALSE)
        }

        preferred_cols <- if (identical(kind, "m/z")) {
            c("mz", "mz.data.A", "mz.data.B", "theoretical_mz", "mz_isotope", "mz_annotated", "mz_annotated_isotope")
        } else {
            c("time", "rt", "retention_time", "time.data.A", "time.data.B", "time_sample", "time_sample_isotope")
        }
        preferred_cols <- intersect(preferred_cols, available_cols)
        candidate_cols <- unique(c(preferred_cols, available_cols))

        best_col <- NULL
        best_values <- NULL
        best_score <- -1L

        for (col_name in candidate_cols) {
            numeric_values <- suppressWarnings(as.numeric(dt[[col_name]]))
            score <- sum(is.finite(numeric_values))
            if (score > best_score) {
                best_score <- score
                best_col <- col_name
                best_values <- numeric_values
            }
        }

        if (is.null(best_col) || best_score <= 0) {
            stop(label, " has non-numeric or missing ", kind, " values in all candidate columns.", call. = FALSE)
        }

        list(name = best_col, values = best_values)
    }

    prepare_match_dt <- function(x, label, mz_name, time_name = NULL) {
        dt <- to_feature_dt(x, label)

        mz_col <- choose_numeric_column(dt, label = label, kind = "m/z")
        dt[, (mz_name) := mz_col$values]
        dt <- dt[is.finite(get(mz_name))]
        if (nrow(dt) == 0) {
            stop(label, " has no valid numeric m/z values after filtering.", call. = FALSE)
        }

        if (!is.null(time_name)) {
            time_col <- choose_numeric_column(
                dt,
                label = label,
                kind = "time",
                excluded = mz_col$name
            )
            dt[, (time_name) := time_col$values]
            dt <- dt[is.finite(get(time_name))]
            if (nrow(dt) == 0) {
                stop(label, " has no valid numeric retention time values after filtering.", call. = FALSE)
            }
        }

        dt
    }

    use_time <- !is.na(time.thresh)
    dataA <- prepare_match_dt(dataA, "dataA", mz_name = "mz.data.A", time_name = if (use_time) "time.data.A" else NULL)
    dataB <- prepare_match_dt(dataB, "dataB", mz_name = "mz.data.B", time_name = if (use_time) "time.data.B" else NULL)

    if (use_time) {
        time.thresh <- suppressWarnings(as.numeric(time.thresh))
        if (!is.finite(time.thresh)) {
            stop("time.thresh must be numeric when provided.", call. = FALSE)
        }
    }

    dataA[, index.A := .original_index__]
    dataB[, index.B := .original_index__]

    dataA[, mz_tol := mz.thresh * mz.data.A / 1000000]
    dataA[, mz_lower := mz.data.A - mz_tol]
    dataA[, mz_upper := mz.data.A + mz_tol]

    data.table::setkey(dataA, mz_lower, mz_upper)
    data.table::setkey(dataB, mz.data.B)

    if (!use_time) {
        matches <- dataA[
            dataB,
            on = list(mz_lower <= mz.data.B, mz_upper >= mz.data.B),
            nomatch = 0L,
            list(index.A = index.A, mz.data.A = mz.data.A, index.B = i.index.B, mz.data.B = i.mz.data.B)
        ]
        return(as.data.frame(matches))
    }

    matches <- dataA[
        dataB,
        on = list(mz_lower <= mz.data.B, mz_upper >= mz.data.B),
        nomatch = 0L,
        list(
            index.A = index.A,
            mz.data.A = mz.data.A,
            time.data.A = time.data.A,
            index.B = i.index.B,
            mz.data.B = i.mz.data.B,
            time.data.B = i.time.data.B
        )
    ]
    matches[, time.difference := abs(time.data.A - time.data.B)]
    final_matches <- matches[time.difference < time.thresh,
        list(index.A, mz.data.A, time.data.A, index.B, mz.data.B, time.data.B, time.difference)
    ]

    as.data.frame(final_matches)
}
