# data-raw/validate_stan_models.R
# Compile Stan models, run on full NE25 validation data, compare to Mplus scores.
#
# Validation criteria:
#   - F scores: correlation > 0.98 with Mplus
#   - Age gradient: positive (older children score higher)
#   - CSEMs: plausible range

library(cmdstanr)

# --- Load data ----------------------------------------------------------------

load("data/item_params.rda")
load("data-raw/validation/validation_ne25_items.rda")
load("data-raw/validation/validation_ne25_ages.rda")
load("data-raw/validation/validation_mplus_scores_f.rda")
load("data-raw/validation/validation_mplus_scores_bifactor.rda")

# Uppercase validation data column names to match Mplus convention in item_params
names(validation_ne25_items) <- toupper(names(validation_ne25_items))

# --- Helper: prepare response matrix ------------------------------------------

prepare_response_matrix <- function(items_df, item_names, n_cat) {
  # Extract item columns and convert to Stan 1-indexed (1, ..., K).
  # Mplus data is 0-indexed but some items in the raw .dat may be offset
  # (e.g., coded 1-4 instead of 0-3). Normalize by subtracting per-item
  # minimum observed value, then add 1. Validate against n_cat.
  y <- as.matrix(items_df[, item_names, drop = FALSE])
  for (j in seq_len(ncol(y))) {
    col_vals <- y[, j]
    col_vals <- col_vals[!is.na(col_vals)]
    if (length(col_vals) == 0) next
    item_min <- min(col_vals)
    if (item_min != 0) {
      y[, j] <- y[, j] - item_min
    }
  }
  # Now 0-indexed. Convert to 1-indexed.
  y <- y + 1L
  # Validate: max response per item should not exceed n_cat
  for (j in seq_len(ncol(y))) {
    col_vals <- y[, j]
    col_vals <- col_vals[!is.na(col_vals)]
    if (length(col_vals) == 0) next
    if (max(col_vals) > n_cat[j]) {
      stop(sprintf(
        "Item %s: max response %d exceeds n_cat %d",
        item_names[j], max(col_vals), n_cat[j]
      ))
    }
  }
  # NA -> 0 (missing indicator for Stan)
  y[is.na(y)] <- 0L
  y
}

# --- Helper: build design matrix ----------------------------------------------

build_design_matrix <- function(ages) {
  cbind(1, log(ages + 1), ages)
}

# =============================================================================
# Part 1: Unidimensional (F factor, 220 items)
# =============================================================================

cat("=== UNIDIMENSIONAL MODEL ===\n\n")

# Identify F-only items (items that load on F but not on any subfactor)
f_items <- item_params$slopes$item[item_params$slopes$factor == "F"]
subfactor_items <- unique(
  item_params$slopes$item[item_params$slopes$factor != "F"]
)
f_only_items <- setdiff(f_items, subfactor_items)
cat("F-only items:", length(f_only_items), "\n")

# Item names are UPPER CASE (from Mplus); validation data columns are also upper
# Verify overlap
available <- intersect(f_only_items, names(validation_ne25_items))
cat("Available in validation data:", length(available), "\n")

# Use only available items
f_only_items <- available

# Get slopes and thresholds for F-only items
f_slopes <- item_params$slopes[
  item_params$slopes$factor == "F" & item_params$slopes$item %in% f_only_items,
]
f_slopes <- f_slopes[match(f_only_items, f_slopes$item), ]

f_thresh <- item_params$thresholds[item_params$thresholds$item %in% f_only_items, ]

# Build n_cat vector
f_n_cat <- item_params$n_cat[f_only_items]

# Build threshold matrix (J x max_cat-1, padded with 0)
f_max_cat <- max(f_n_cat)
J_f <- length(f_only_items)
f_d_matrix <- matrix(0, nrow = J_f, ncol = f_max_cat - 1)
for (idx in seq_along(f_only_items)) {
  item <- f_only_items[idx]
  item_thresh <- f_thresh[f_thresh$item == item, ]
  item_thresh <- item_thresh[order(item_thresh$threshold), ]
  for (row_idx in seq_len(nrow(item_thresh))) {
    f_d_matrix[idx, item_thresh$threshold[row_idx]] <- item_thresh$d[row_idx]
  }
}

# Response matrix
# The F-only model uses 2785 persons but Mplus F scores have 2781 rows
# We'll run on all 2785 and compare the subset
y_f <- prepare_response_matrix(validation_ne25_items, f_only_items, f_n_cat)
X_f <- build_design_matrix(validation_ne25_ages$years_old)

cat("Response matrix:", nrow(y_f), "x", ncol(y_f), "\n")
cat("Max categories:", f_max_cat, "\n")

# --- Compile unidimensional model ---------------------------------------------

cat("\nCompiling unidimensional model...\n")
mod_f <- cmdstan_model(
  "inst/stan/grm_unidimensional.stan",
)
cat("Compilation successful.\n")

# --- Prepare Stan data --------------------------------------------------------

stan_data_f <- list(
  N = nrow(y_f),
  J = ncol(y_f),
  max_cat = f_max_cat,
  n_cat = as.array(as.integer(f_n_cat)),
  y = y_f,
  a = f_slopes$slope,
  d = f_d_matrix,
  X = X_f
)

# --- Run BFGS optimization ---------------------------------------------------

cat("\nRunning BFGS optimization (unidimensional)...\n")
fit_f <- mod_f$optimize(
  data = stan_data_f,
  algorithm = "lbfgs",
  history_size = 50,
  init = 0.1,
  iter = 20000,
  tol_rel_grad = 1e-8,
  refresh = 500
)

cat("Optimization complete.\n")

# --- Extract theta estimates --------------------------------------------------

theta_names <- paste0("theta[", 1:stan_data_f$N, "]")
theta_f <- fit_f$summary(variables = "theta")
cat("\nTheta summary (first 5):\n")
print(head(theta_f, 5))

# --- Compare to Mplus F scores ------------------------------------------------

# Mplus F scores have 2781 rows; our model has 2785
# The Mplus kidsights model dropped 4 cases; bifactor model kept all 2785
# Use bifactor F scores (2785 rows) for comparison since row counts match
our_f_scores <- theta_f$estimate
mplus_f_scores <- validation_mplus_scores_bifactor$F

cat("\n=== VALIDATION: F SCORES ===\n")
cat(
  "Our scores: N =", length(our_f_scores),
  "range [", round(min(our_f_scores), 3), ",",
  round(max(our_f_scores), 3), "]\n"
)
cat(
  "Mplus scores: N =", length(mplus_f_scores),
  "range [", round(min(mplus_f_scores), 3), ",",
  round(max(mplus_f_scores), 3), "]\n"
)

r_f <- cor(our_f_scores, mplus_f_scores)
cat("Correlation:", round(r_f, 5), "\n")
cat("Correlation > 0.98:", r_f > 0.98, "\n")

rmse_f <- sqrt(mean((our_f_scores - mplus_f_scores)^2))
cat("RMSE:", round(rmse_f, 4), "\n")

# Age gradient
age_cor <- cor(validation_ne25_ages$years_old, our_f_scores)
cat(
  "Age-score correlation:", round(age_cor, 4),
  "(positive =", age_cor > 0, ")\n"
)

# Regression coefficients
beta_f <- fit_f$summary(variables = "beta")
cat("\nRegression coefficients (beta):\n")
print(beta_f)

# =============================================================================
# Part 2: Bifactor (GEN + 5 subfactors, 44 PS items)
# =============================================================================

cat("\n\n=== BIFACTOR MODEL ===\n\n")

# Identify PS items and their factor loadings
ps_factors <- c("EAT", "EXT", "INT", "SLE", "SOC")
gen_slopes <- item_params$slopes[item_params$slopes$factor == "GEN", ]
sub_slopes <- item_params$slopes[item_params$slopes$factor %in% ps_factors, ]

# PS items = items that appear in GEN loadings
ps_items <- gen_slopes$item
cat("PS items:", length(ps_items), "\n")

# Verify all PS items also have a subfactor loading
ps_in_sub <- ps_items %in% sub_slopes$item
cat("PS items with subfactor loading:", sum(ps_in_sub), "/", length(ps_items), "\n")

# Available in validation data
ps_available <- intersect(ps_items, names(validation_ne25_items))
cat("Available in validation data:", length(ps_available), "\n")
ps_items <- ps_available

# Order consistently
gen_slopes <- gen_slopes[match(ps_items, gen_slopes$item), ]

# Build subfactor index and slopes
# For each PS item, find which subfactor it loads on
sub_idx <- integer(length(ps_items))
a_sub_vec <- numeric(length(ps_items))
for (idx in seq_along(ps_items)) {
  item <- ps_items[idx]
  row <- sub_slopes[sub_slopes$item == item, ]
  sub_idx[idx] <- match(row$factor, ps_factors)
  a_sub_vec[idx] <- row$slope
}

cat("Subfactor distribution:\n")
print(table(ps_factors[sub_idx]))

# Thresholds
ps_thresh <- item_params$thresholds[item_params$thresholds$item %in% ps_items, ]
ps_n_cat <- item_params$n_cat[ps_items]
ps_max_cat <- max(ps_n_cat)
J_ps <- length(ps_items)

ps_d_matrix <- matrix(0, nrow = J_ps, ncol = ps_max_cat - 1)
for (idx in seq_along(ps_items)) {
  item <- ps_items[idx]
  item_thresh <- ps_thresh[ps_thresh$item == item, ]
  item_thresh <- item_thresh[order(item_thresh$threshold), ]
  for (row_idx in seq_len(nrow(item_thresh))) {
    ps_d_matrix[idx, item_thresh$threshold[row_idx]] <- item_thresh$d[row_idx]
  }
}

# Response matrix
y_ps <- prepare_response_matrix(validation_ne25_items, ps_items, ps_n_cat)
X_ps <- build_design_matrix(validation_ne25_ages$years_old)

cat("Response matrix:", nrow(y_ps), "x", ncol(y_ps), "\n")
cat("Max categories:", ps_max_cat, "\n")

# --- Compile bifactor model ---------------------------------------------------

cat("\nCompiling bifactor model...\n")
mod_bf <- cmdstan_model(
  "inst/stan/grm_bifactor.stan",
)
cat("Compilation successful.\n")

# --- Prepare Stan data --------------------------------------------------------

stan_data_bf <- list(
  N = nrow(y_ps),
  J = ncol(y_ps),
  K_sub = length(ps_factors),
  max_cat = ps_max_cat,
  n_cat = as.array(as.integer(ps_n_cat)),
  y = y_ps,
  a_gen = gen_slopes$slope,
  a_sub = a_sub_vec,
  sub_idx = as.array(sub_idx),
  d = ps_d_matrix,
  X = X_ps
)

# --- Run BFGS optimization ---------------------------------------------------

cat("\nRunning BFGS optimization (bifactor)...\n")
fit_bf <- mod_bf$optimize(
  data = stan_data_bf,
  algorithm = "lbfgs",
  history_size = 50,
  init = 0.1,
  iter = 20000,
  tol_rel_grad = 1e-8,
  refresh = 500
)

cat("Optimization complete.\n")

# --- Extract and compare GEN scores ------------------------------------------

theta_gen <- fit_bf$summary(variables = "theta_gen")
our_gen_scores <- theta_gen$estimate

mplus_gen_scores <- validation_mplus_scores_bifactor$GEN

cat("\n=== VALIDATION: GEN SCORES ===\n")
cat(
  "Our scores: N =", length(our_gen_scores),
  "range [", round(min(our_gen_scores), 3), ",",
  round(max(our_gen_scores), 3), "]\n"
)
cat(
  "Mplus scores: N =", length(mplus_gen_scores),
  "range [", round(min(mplus_gen_scores), 3), ",",
  round(max(mplus_gen_scores), 3), "]\n"
)

r_gen <- cor(our_gen_scores, mplus_gen_scores)
cat("Correlation:", round(r_gen, 5), "\n")

# --- Compare subfactor scores ------------------------------------------------

cat("\n=== VALIDATION: SUBFACTOR SCORES ===\n")
for (k in seq_along(ps_factors)) {
  fac <- ps_factors[k]
  theta_sub_k <- fit_bf$summary(
    variables = paste0("theta_sub[", 1:stan_data_bf$N, ",", k, "]")
  )
  our_sub <- theta_sub_k$estimate
  mplus_sub <- validation_mplus_scores_bifactor[[fac]]

  r_sub <- cor(our_sub, mplus_sub)
  cat(sprintf(
    "  %s: r = %.5f, our range [%.3f, %.3f], Mplus range [%.3f, %.3f]\n",
    fac, r_sub,
    min(our_sub), max(our_sub),
    min(mplus_sub), max(mplus_sub)
  ))
}

# --- Final summary ------------------------------------------------------------

cat("\n\n=== VALIDATION SUMMARY ===\n")
cat(sprintf("F scores:   r = %.5f  (target > 0.98: %s)\n", r_f, r_f > 0.98))
cat(sprintf("GEN scores: r = %.5f\n", r_gen))
cat(sprintf("Age-F cor:  %.4f  (positive: %s)\n", age_cor, age_cor > 0))
cat(sprintf("F RMSE:     %.4f\n", rmse_f))
