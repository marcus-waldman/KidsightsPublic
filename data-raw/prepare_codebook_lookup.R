# data-raw/prepare_codebook_lookup.R
# Flatten codebook.json into a long lexicon lookup table.
#
# Each row maps: equate_id (canonical) x lexicon -> column_name
# Used by detect_lexicon() to auto-identify which lexicon a dataset uses.

library(jsonlite)

# --- Configuration -----------------------------------------------------------

codebook_path <- file.path(
  "C:/Users/marcu/git-repositories/Kidsights-Data-Platform",
  "codebook/data/codebook.json"
)

# --- Read codebook ------------------------------------------------------------

cb <- fromJSON(codebook_path)
items <- cb$items

# --- Flatten into long table --------------------------------------------------

rows <- vector("list", length(items) * 10)
idx <- 0L

for (equate_id in names(items)) {
  lexicons <- items[[equate_id]]$lexicons
  if (is.null(lexicons)) next
  for (lex_name in names(lexicons)) {
    col_name <- lexicons[[lex_name]]
    if (is.null(col_name) || length(col_name) == 0) next
    if (any(is.na(col_name))) next
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      equate_id = equate_id,
      lexicon = lex_name,
      column_name = col_name,
      stringsAsFactors = FALSE
    )
  }
}

codebook_lookup <- do.call(rbind, rows[seq_len(idx)])
rownames(codebook_lookup) <- NULL

# --- Save ---------------------------------------------------------------------

usethis::use_data(codebook_lookup, overwrite = TRUE)

cat("Saved codebook_lookup.rda\n")
cat("  Rows:", nrow(codebook_lookup), "\n")
cat("  Unique items:", length(unique(codebook_lookup$equate_id)), "\n")
cat("  Lexicons:", paste(sort(unique(codebook_lookup$lexicon)), collapse = ", "), "\n")
