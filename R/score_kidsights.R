#' Score Kidsights (F factor)
#'
#' Compute MAP estimates of the Kidsights developmental score using a
#' Graded Response Model with fixed item parameters and an age-informed prior.
#' Uses L-BFGS optimization via CmdStan.
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
#' @return A data frame with the ID columns and a \code{theta} column
#'   containing MAP score estimates. Observations with fewer than
#'   \code{min_responses} items have \code{theta = NA}.
#' @export
score_kidsights <- function(data, id_cols, min_responses = 5, ...) {
  validate_inputs(data, id_cols)

  lex <- detect_lexicon(data)
  mapping <- lex$mapping

  item_params <- KidsightsPublic::item_params
  f_items <- item_params$slopes$item[item_params$slopes$factor == "F"]
  subfactor_items <- unique(
    item_params$slopes$item[item_params$slopes$factor != "F"]
  )
  f_only_items <- setdiff(f_items, subfactor_items)

  mapped_equate <- toupper(mapping)
  available_mask <- mapped_equate %in% f_only_items
  mapping <- mapping[available_mask]

  if (length(mapping) == 0) {
    stop(
      "No F-factor items found in the data after lexicon mapping.",
      call. = FALSE
    )
  }

  item_data <- data[, names(mapping), drop = FALSE]
  names(item_data) <- toupper(mapping)
  used_items <- toupper(mapping)

  item_data <- validate_responses(item_data, used_items, item_params$n_cat)

  y <- as.matrix(item_data)
  n_responses <- rowSums(!is.na(y))
  scoreable <- n_responses >= min_responses
  y_stan <- y[scoreable, , drop = FALSE] + 1L
  y_stan[is.na(y_stan)] <- 0L

  ages <- data$years_old[scoreable]
  x_mat <- cbind(1, log(ages + 1), ages)

  slopes <- item_params$slopes[
    item_params$slopes$factor == "F" & item_params$slopes$item %in% used_items,
  ]
  slopes <- slopes[match(used_items, slopes$item), ]

  n_cat_vec <- item_params$n_cat[used_items]
  max_cat <- max(n_cat_vec)
  j_items <- length(used_items)

  d_matrix <- build_threshold_matrix(item_params$thresholds, used_items,
                                     n_cat_vec, max_cat)

  stan_data <- list(
    N = nrow(y_stan),
    J = j_items,
    max_cat = max_cat,
    n_cat = as.array(as.integer(n_cat_vec)),
    y = y_stan,
    a = slopes$slope,
    d = d_matrix,
    X = x_mat
  )

  stan_file <- system.file("stan", "grm_unidimensional.stan",
                           package = "KidsightsPublic")
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

  theta_all <- rep(NA_real_, nrow(data))
  theta_est <- fit$summary(variables = "theta")$estimate
  theta_all[scoreable] <- theta_est

  result <- data[, id_cols, drop = FALSE]
  result$theta <- theta_all
  result
}
