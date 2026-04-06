# data-raw/prepare_validation_data.R
# Pull NE25 validation data from Kidsights-Data-Platform.
# Saves to data-raw/validation/ (excluded from package build via .Rbuildignore).
#
# Uses MplusAutomation::readModels() to read SAVEDATA with proper column names.

library(MplusAutomation)

# --- Configuration -----------------------------------------------------------

mplus_dir <- file.path(
  "C:/Users/marcu/git-repositories/Kidsights-Data-Platform",
  "calibration/ne25/manual_2023_scale/mplus"
)

output_dir <- "data-raw/validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Parse column names from mplus_dat.inp -----------------------------------

inp_lines <- readLines(file.path(mplus_dir, "mplus_dat.inp"))
names_text <- paste(inp_lines, collapse = " ")
names_text <- sub(".*NAMES\\s*=\\s*", "", names_text, ignore.case = TRUE)
names_text <- sub("\\s*;.*", "", names_text)
col_names <- unlist(strsplit(trimws(names_text), "\\s+"))

# --- Load item responses from mplus_dat.dat -----------------------------------

dat <- read.table(
  file.path(mplus_dir, "mplus_dat.dat"),
  col.names = col_names,
  na.strings = "."
)
cat("Loaded mplus_dat.dat:", nrow(dat), "rows x", ncol(dat), "cols\n")

# Identify item columns (everything after the 14 id/covariate columns)
id_cov_names <- c(
  "rid", "pid", "recordid",
  "logyrs", "yrs3", "school", "logfpl", "phq2",
  "schXyrs3", "fplXyrs3", "phqXyrs3", "black", "hisp", "other"
)
item_names <- setdiff(col_names, id_cov_names)

# Extract item responses
validation_ne25_items <- dat[, item_names]
cat("Item matrix:", nrow(validation_ne25_items), "x", ncol(validation_ne25_items), "\n")

# --- Back-compute age from logyrs ---------------------------------------------
# logyrs = ln(years_old + 1), so years_old = exp(logyrs) - 1

validation_ne25_ages <- data.frame(
  rid = dat$rid,
  years_old = exp(dat$logyrs) - 1
)
cat("Age range:", round(min(validation_ne25_ages$years_old), 2), "to",
    round(max(validation_ne25_ages$years_old), 2), "years\n")

# --- Load Mplus ground truth scores via MplusAutomation -----------------------

# Bifactor model (all 7 factors + F)
mod_all <- readModels(file.path(mplus_dir, "all_2023_calibration_ne25.out"))
savedata_all <- mod_all$savedata

# Factor score columns: F, F_SE, EAT, EAT_SE, ..., GEN, GEN_SE, RID
factor_names <- c("F", "EAT", "EXT", "INT", "SLE", "SOC", "GEN")
theta_cols <- factor_names
se_cols <- paste0(factor_names, "_SE")

validation_mplus_scores_bifactor <- savedata_all[, c(theta_cols, se_cols)]
cat("Bifactor scores:", nrow(validation_mplus_scores_bifactor), "rows,",
    ncol(validation_mplus_scores_bifactor), "cols\n")

# Kidsights (F only) model
mod_k <- readModels(file.path(mplus_dir, "kidsights_2023_calibration_ne25.out"))
savedata_k <- mod_k$savedata

validation_mplus_scores_f <- savedata_k[, c("F", "F_SE")]
cat("F scores:", nrow(validation_mplus_scores_f), "rows\n")
cat("  Note: F-only model has", nrow(savedata_k), "rows vs bifactor",
    nrow(savedata_all), "rows\n")

# --- Save ---------------------------------------------------------------------

save(validation_ne25_items,
     file = file.path(output_dir, "validation_ne25_items.rda"),
     compress = "xz")
save(validation_ne25_ages,
     file = file.path(output_dir, "validation_ne25_ages.rda"),
     compress = "xz")
save(validation_mplus_scores_f,
     file = file.path(output_dir, "validation_mplus_scores_f.rda"),
     compress = "xz")
save(validation_mplus_scores_bifactor,
     file = file.path(output_dir, "validation_mplus_scores_bifactor.rda"),
     compress = "xz")

cat("\nSaved validation data to", output_dir, "\n")
