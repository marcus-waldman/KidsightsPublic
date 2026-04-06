#' Score Psychosocial Domains (Bifactor Model)
#'
#' Compute MAP estimates of psychosocial domain scores using a bifactor
#' Graded Response Model with fixed item parameters and per-factor
#' age-informed priors. Returns a general psychosocial factor (GEN) plus
#' five subfactor scores (EAT, EXT, INT, SLE, SOC).
#'
#' @param data A data frame containing item responses, a \code{years_old}
#'   column with age in years, and one or more ID columns.
#' @param id_cols Character vector of column names that uniquely identify
#'   each observation.
#' @param min_responses Minimum number of non-missing item responses required
#'   to score an observation. Observations with fewer responses receive
#'   \code{NA}. Defaults to 5.
#' @param ... Additional arguments passed to \code{cmdstanr} model
#'   \code{$optimize()}.
#' @return A data frame with the ID columns and score columns:
#'   \code{theta_gen}, \code{theta_eat}, \code{theta_ext}, \code{theta_int},
#'   \code{theta_sle}, \code{theta_soc}.
#'   Observations with fewer than \code{min_responses} items have \code{NA}.
#' @export
score_psychosocial <- function(data, id_cols, min_responses = 5, ...) {
  validate_inputs(data, id_cols)

  lex <- detect_lexicon(data)
  mapping <- lex$mapping

  item_params <- KidsightsPublic::item_params
  ps_factors <- c("EAT", "EXT", "INT", "SLE", "SOC")

  gen_slopes <- item_params$slopes[item_params$slopes$factor == "GEN", ]
  sub_slopes <- item_params$slopes[item_params$slopes$factor %in% ps_factors, ]
  ps_items <- gen_slopes$item

  mapped_equate <- toupper(mapping)
  available_mask <- mapped_equate %in% ps_items
  mapping <- mapping[available_mask]

  if (length(mapping) == 0) {
    stop(
      "No psychosocial items found in the data after lexicon mapping.",
      call. = FALSE
    )
  }

  item_data <- data[, names(mapping), drop = FALSE]
  names(item_data) <- toupper(mapping)
  used_items <- toupper(mapping)

  item_data <- validate_responses(item_data, used_items, item_params$n_cat)

  gen_slopes <- gen_slopes[gen_slopes$item %in% used_items, ]
  gen_slopes <- gen_slopes[match(used_items, gen_slopes$item), ]

  sub_idx <- integer(length(used_items))
  a_sub_vec <- numeric(length(used_items))
  for (idx in seq_along(used_items)) {
    row <- sub_slopes[sub_slopes$item == used_items[idx], ]
    sub_idx[idx] <- match(row$factor, ps_factors)
    a_sub_vec[idx] <- row$slope
  }

  y <- as.matrix(item_data)
  n_responses <- rowSums(!is.na(y))
  scoreable <- n_responses >= min_responses
  y_stan <- y[scoreable, , drop = FALSE] + 1L
  y_stan[is.na(y_stan)] <- 0L

  ages <- data$years_old[scoreable]
  x_mat <- cbind(1, log(ages + 1), ages)

  n_cat_vec <- item_params$n_cat[used_items]
  max_cat <- max(n_cat_vec)

  d_matrix <- build_threshold_matrix(
    item_params$thresholds, used_items,
    n_cat_vec, max_cat
  )

  stan_data <- list(
    N = nrow(y_stan),
    J = length(used_items),
    K_sub = length(ps_factors),
    max_cat = max_cat,
    n_cat = as.array(as.integer(n_cat_vec)),
    y = y_stan,
    a_gen = gen_slopes$slope,
    a_sub = a_sub_vec,
    sub_idx = as.array(sub_idx),
    d = d_matrix,
    X = x_mat
  )

  stan_file <- system.file("stan", "grm_bifactor.stan",
    package = "KidsightsPublic"
  )
  mod <- cmdstanr::cmdstan_model(stan_file)

  fit <- mod$optimize(
    data = stan_data,
    algorithm = "lbfgs",
    history_size = 50,
    init = 0.1,
    iter = 20000,
    tol_rel_grad = 1e-8,
    refresh = 0,
    ...
  )

  n_scoreable <- sum(scoreable)
  result <- data[, id_cols, drop = FALSE]

  theta_gen_est <- fit$summary(variables = "theta_gen")$estimate
  result$theta_gen <- NA_real_
  result$theta_gen[scoreable] <- theta_gen_est

  for (k in seq_along(ps_factors)) {
    var_names <- paste0("theta_sub[", seq_len(n_scoreable), ",", k, "]")
    theta_sub_est <- fit$summary(variables = var_names)$estimate
    col_name <- paste0("theta_", tolower(ps_factors[k]))
    result[[col_name]] <- NA_real_
    result[[col_name]][scoreable] <- theta_sub_est
  }

  result
}
