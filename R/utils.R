#' Validate common inputs for scoring functions
#' @param data Data frame.
#' @param id_cols Character vector of ID column names.
#' @noRd
validate_inputs <- function(data, id_cols) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  missing_ids <- setdiff(id_cols, names(data))
  if (length(missing_ids) > 0) {
    stop(
      "ID columns not found in data: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  id_data <- data[, id_cols, drop = FALSE]
  if (anyDuplicated(id_data)) {
    stop(
      "ID columns do not uniquely identify observations. ",
      "Duplicate rows found for: ",
      paste(id_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!("years_old" %in% names(data))) {
    stop(
      "Column `years_old` not found in data. ",
      "Age must be provided in raw years.",
      call. = FALSE
    )
  }

  if (any(is.na(data$years_old))) {
    stop("`years_old` contains missing values.", call. = FALSE)
  }

  if (any(data$years_old < 0)) {
    stop("`years_old` contains negative values.", call. = FALSE)
  }
}


#' Validate item responses against expected categories
#'
#' Responses outside \code{[0, n_cat - 1]} are set to \code{NA} with a
#' warning. This handles items whose raw coding is shifted relative to the
#' calibration (e.g., coded 1-4 instead of 0-3).
#'
#' @param item_data Data frame of item responses (columns named as equate IDs).
#'   Modified in place (caller should use the returned value).
#' @param item_names Character vector of item names.
#' @param n_cat Named integer vector of category counts per item.
#' @return The \code{item_data} data frame with out-of-range values set to
#'   \code{NA}.
#' @noRd
validate_responses <- function(item_data, item_names, n_cat) {
  flagged_items <- character(0)
  for (item in item_names) {
    vals <- item_data[[item]]
    out_of_range <- !is.na(vals) & (vals < 0 | vals > n_cat[item] - 1L)
    if (any(out_of_range)) {
      flagged_items <- c(flagged_items, item)
      item_data[[item]][out_of_range] <- NA
    }
  }
  if (length(flagged_items) > 0) {
    warning(
      "Responses outside valid range set to NA for items: ",
      paste(flagged_items, collapse = ", "),
      call. = FALSE
    )
  }
  item_data
}


#' Build threshold matrix for Stan
#' @param thresholds Data frame with columns: item, threshold, d.
#' @param item_names Character vector of item names (defines row order).
#' @param n_cat Named integer vector of category counts.
#' @param max_cat Maximum number of categories.
#' @return Matrix of dimension \code{length(item_names)} by \code{max_cat - 1},
#'   padded with zeros.
#' @noRd
build_threshold_matrix <- function(thresholds, item_names, n_cat, max_cat) {
  j <- length(item_names)
  d_matrix <- matrix(0, nrow = j, ncol = max_cat - 1)
  for (idx in seq_along(item_names)) {
    item <- item_names[idx]
    item_thresh <- thresholds[thresholds$item == item, ]
    item_thresh <- item_thresh[order(item_thresh$threshold), ]
    for (row_idx in seq_len(nrow(item_thresh))) {
      d_matrix[idx, item_thresh$threshold[row_idx]] <- item_thresh$d[row_idx]
    }
  }
  d_matrix
}
