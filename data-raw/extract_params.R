# data-raw/extract_params.R
# Extract item parameters from Mplus calibration output using MplusAutomation.
#
# Source: all_2023_calibration_ne25.out (bifactor model)
# Contains all 7 factors: F (Kidsights overall), GEN, EAT, EXT, INT, SLE, SOC
#
# CRITICAL: Mplus threshold tau is sign-flipped to Stan intercept d = -tau.

library(MplusAutomation)

# --- Configuration -----------------------------------------------------------

mplus_out_path <- file.path(
  "C:/Users/marcu/git-repositories/Kidsights-Data-Platform",
  "calibration/ne25/manual_2023_scale/mplus",
  "all_2023_calibration_ne25.out"
)

# --- Read Mplus output --------------------------------------------------------

mod <- readModels(mplus_out_path)
unstd <- mod$parameters$unstandardized

# --- Extract slopes (BY rows) ------------------------------------------------

slopes_raw <- unstd[grepl("\\.BY$", unstd$paramHeader), ]
slopes <- data.frame(
  item = slopes_raw$param,
  factor = sub("\\.BY$", "", slopes_raw$paramHeader),
  slope = slopes_raw$est,
  stringsAsFactors = FALSE
)

# --- Extract thresholds (sign-flip: d = -tau) ---------------------------------

thresh_raw <- unstd[unstd$paramHeader == "Thresholds", ]

# Parse item name and threshold index from param (e.g., "AA4$1" -> "AA4", 1)
thresh_parts <- strsplit(thresh_raw$param, "\\$")
thresh <- data.frame(
  item = vapply(thresh_parts, `[`, character(1), 1),
  threshold = as.integer(vapply(thresh_parts, `[`, character(1), 2)),
  d = -thresh_raw$est,
  stringsAsFactors = FALSE
)

# --- Compute number of categories per item ------------------------------------

# n_cat = max threshold index + 1 (e.g., 1 threshold -> 2 categories)
max_thresh <- tapply(thresh$threshold, thresh$item, max)
n_cat <- max_thresh + 1L
n_cat <- n_cat[order(names(n_cat))]

# --- Factor metadata ----------------------------------------------------------

factors <- sort(unique(slopes$factor))

# --- Assemble and save --------------------------------------------------------

item_params <- list(
  slopes = slopes,
  thresholds = thresh,
  n_cat = n_cat,
  factors = factors
)

usethis::use_data(item_params, overwrite = TRUE)

cat("Saved item_params.rda\n")
cat("  Slopes:", nrow(slopes), "rows across", length(factors), "factors\n")
cat("  Thresholds:", nrow(thresh), "rows (sign-flipped: d = -tau)\n")
cat("  Items:", length(n_cat), "\n")
cat("  Factors:", paste(factors, collapse = ", "), "\n")
