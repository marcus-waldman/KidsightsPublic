#' Detect which lexicon a dataset uses
#'
#' Checks column names of a data frame against all lexicons in the bundled
#' codebook lookup table. Returns the best-matching lexicon and a mapping
#' from data column names to canonical equate IDs.
#'
#' @param data A data frame with item response columns.
#' @param min_overlap Minimum fraction of lexicon items that must be present
#'   in the data for a match. Defaults to 0.1.
#' @return A list with components:
#'   \describe{
#'     \item{lexicon}{Character string: name of the best-matching lexicon.}
#'     \item{mapping}{Named character vector: names are data column names,
#'       values are canonical equate IDs.}
#'     \item{n_matched}{Integer: number of items matched.}
#'     \item{n_lexicon}{Integer: total items in the matched lexicon.}
#'   }
#' @export
detect_lexicon <- function(data, min_overlap = 0.1) {
  lookup <- KidsightsPublic::codebook_lookup

  data_cols_upper <- toupper(names(data))
  lexicons <- unique(lookup$lexicon)

  best_lexicon <- NULL
  best_count <- 0L
  best_n_total <- 0L

  for (lex in lexicons) {
    lex_rows <- lookup[lookup$lexicon == lex, ]
    lex_cols_upper <- toupper(lex_rows$column_name)
    matched <- sum(data_cols_upper %in% lex_cols_upper)
    if (matched > best_count) {
      best_count <- matched
      best_lexicon <- lex
      best_n_total <- nrow(lex_rows)
    }
  }

  if (is.null(best_lexicon) || best_count / best_n_total < min_overlap) {
    stop(
      "No lexicon matched the data columns. ",
      "Best match: ", if (is.null(best_lexicon)) "none" else best_lexicon,
      " (", best_count, "/", best_n_total, " items). ",
      "Ensure item column names follow a recognized lexicon.",
      call. = FALSE
    )
  }

  lex_rows <- lookup[lookup$lexicon == best_lexicon, ]
  lex_cols_upper <- toupper(lex_rows$column_name)

  matched_idx <- which(data_cols_upper %in% lex_cols_upper)
  mapping <- character(length(matched_idx))
  names(mapping) <- names(data)[matched_idx]

  for (i in seq_along(matched_idx)) {
    data_col_upper <- data_cols_upper[matched_idx[i]]
    row_idx <- which(lex_cols_upper == data_col_upper)[1]
    mapping[i] <- lex_rows$equate_id[row_idx]
  }

  list(
    lexicon = best_lexicon,
    mapping = mapping,
    n_matched = best_count,
    n_lexicon = best_n_total
  )
}
