#' Simulate complete GRM item responses
#'
#' Generate simulated item responses from the Graded Response Model using
#' the package's fixed item parameters. Produces complete (no missingness)
#' 0-indexed responses for all 264 Kidsights items (220 F-factor plus
#' 44 psychosocial).
#'
#' When \code{theta} is \code{NULL}, latent trait values are drawn from an
#' age-informed prior estimated from the NE25 validation sample (N = 2,785).
#' The prior uses the same design matrix as the scoring models:
#' \code{theta ~ N(b0 + b1 * ln(age + 1) + b2 * age, sigma)}.
#'
#' @param years_old Numeric scalar or vector of ages in years. If scalar and
#'   \code{n > 1}, the value is recycled. If a vector, its length determines
#'   the number of observations.
#' @param n Integer number of observations to simulate. Used only when
#'   \code{years_old} is a scalar. Defaults to \code{1L}.
#' @param theta Optional numeric matrix of latent trait values with
#'   \code{n} rows and 7 named columns:
#'   \code{c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")}.
#'   If \code{NULL} (the default), thetas are drawn from the age-informed
#'   prior.
#' @param seed Optional integer seed for reproducibility. The random number
#'   generator state is saved on entry and restored on exit so the caller's
#'   RNG stream is not affected.
#' @return A data frame with \code{n} rows and columns:
#'   \describe{
#'     \item{id}{Integer observation identifier (1 to n).}
#'     \item{years_old}{Age in years.}
#'     \item{...}{One column per item (264 columns), named by canonical
#'       equate ID, containing 0-indexed integer responses.}
#'   }
#' @export
#' @examples
#' # Simulate 5 observations at age 3
#' sim <- simulate_responses(years_old = 3, n = 5, seed = 42)
#' head(sim[, 1:6])
#'
#' # Simulate with varying ages
#' sim2 <- simulate_responses(years_old = c(1, 2, 3, 4, 5), seed = 123)
#'
#' \donttest{
#' # Supply custom theta values
#' theta_mat <- matrix(0,
#'   nrow = 5, ncol = 7,
#'   dimnames = list(
#'     NULL,
#'     c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")
#'   )
#' )
#' sim3 <- simulate_responses(years_old = 3, n = 5, theta = theta_mat)
#' }
simulate_responses <- function(years_old, n = 1L, theta = NULL, seed = NULL) {
  # --- Input validation -------------------------------------------------------
  if (!is.numeric(years_old)) {
    stop("`years_old` must be numeric.", call. = FALSE)
  }
  if (any(is.na(years_old))) {
    stop("`years_old` contains missing values.", call. = FALSE)
  }
  if (any(years_old < 0)) {
    stop("`years_old` contains negative values.", call. = FALSE)
  }

  if (length(years_old) > 1) {
    n <- length(years_old)
  } else {
    n <- as.integer(n)
    years_old <- rep(years_old, n)
  }

  # --- Seed handling ----------------------------------------------------------
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit(
      {
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = globalenv())
        } else {
          assign(".Random.seed", old_seed, envir = globalenv())
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }

  # --- Load item parameters ---------------------------------------------------
  ip <- KidsightsPublic::item_params

  # F-only items (same logic as score_kidsights)
  f_items <- ip$slopes$item[ip$slopes$factor == "F"]
  subfactor_items <- unique(ip$slopes$item[ip$slopes$factor != "F"])
  f_only_items <- setdiff(f_items, subfactor_items)

  # PS items (same logic as score_psychosocial)
  ps_factors <- c("EAT", "EXT", "INT", "SLE", "SOC")
  gen_slopes <- ip$slopes[ip$slopes$factor == "GEN", ]
  sub_slopes <- ip$slopes[ip$slopes$factor %in% ps_factors, ]
  ps_items <- gen_slopes$item

  all_items <- c(f_only_items, ps_items)

  # --- Theta generation / validation ------------------------------------------
  factor_names <- c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")

  if (is.null(theta)) {
    theta <- generate_theta_from_prior(years_old, n, factor_names)
  } else {
    if (!is.matrix(theta) || !is.numeric(theta)) {
      stop("`theta` must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(theta) != n) {
      stop(
        "`theta` must have ", n, " rows (one per observation).",
        call. = FALSE
      )
    }
    if (ncol(theta) != 7) {
      stop("`theta` must have 7 columns.", call. = FALSE)
    }
    if (is.null(colnames(theta)) ||
      !all(factor_names %in% colnames(theta))) {
      stop(
        "`theta` columns must be named: ",
        paste(factor_names, collapse = ", "),
        call. = FALSE
      )
    }
    theta <- theta[, factor_names, drop = FALSE]
  }

  # --- Build item parameter structures ----------------------------------------
  # F-only slopes (ordered to match f_only_items)
  f_slopes_df <- ip$slopes[
    ip$slopes$factor == "F" & ip$slopes$item %in% f_only_items,
  ]
  f_slopes_df <- f_slopes_df[match(f_only_items, f_slopes_df$item), ]

  # PS slopes (ordered to match ps_items)
  gen_slopes <- gen_slopes[match(ps_items, gen_slopes$item), ]
  ps_sub_idx <- character(length(ps_items))
  ps_a_sub <- numeric(length(ps_items))
  for (idx in seq_along(ps_items)) {
    rows <- sub_slopes[sub_slopes$item == ps_items[idx], ]
    best <- which.max(abs(rows$slope))
    ps_sub_idx[idx] <- rows$factor[best]
    ps_a_sub[idx] <- rows$slope[best]
  }

  # Threshold matrices
  n_cat_f <- ip$n_cat[f_only_items]
  max_cat_f <- max(n_cat_f)
  d_matrix_f <- build_threshold_matrix(
    ip$thresholds, f_only_items, n_cat_f, max_cat_f
  )

  n_cat_ps <- ip$n_cat[ps_items]
  max_cat_ps <- max(n_cat_ps)
  d_matrix_ps <- build_threshold_matrix(
    ip$thresholds, ps_items, n_cat_ps, max_cat_ps
  )

  # --- Simulate F-only item responses -----------------------------------------
  result <- matrix(0L, nrow = n, ncol = length(all_items))

  for (j in seq_along(f_only_items)) {
    eta <- f_slopes_df$slope[j] * theta[, "F"]
    d_vec <- d_matrix_f[j, seq_len(n_cat_f[j] - 1L)]
    result[, j] <- simulate_grm_item(eta, d_vec, n_cat_f[j])
  }

  # --- Simulate PS item responses ---------------------------------------------
  offset <- length(f_only_items)
  for (j in seq_along(ps_items)) {
    eta <- gen_slopes$slope[j] * theta[, "GEN"] +
      ps_a_sub[j] * theta[, ps_sub_idx[j]]
    d_vec <- d_matrix_ps[j, seq_len(n_cat_ps[j] - 1L)]
    result[, offset + j] <- simulate_grm_item(eta, d_vec, n_cat_ps[j])
  }

  # --- Assemble output --------------------------------------------------------
  df <- data.frame(id = seq_len(n), years_old = years_old)
  item_df <- as.data.frame(result)
  names(item_df) <- all_items
  cbind(df, item_df)
}


#' Simulate a single GRM item for n persons
#'
#' @param eta Numeric vector of linear predictor values (length n).
#' @param d_vec Numeric vector of thresholds (length n_cat - 1).
#' @param n_cat Integer number of response categories.
#' @return Integer vector of 0-indexed responses (length n).
#' @noRd
simulate_grm_item <- function(eta, d_vec, n_cat) {
  n <- length(eta)
  # Cumulative probs: P(Y >= k+1) = plogis(eta + d[k]) for k = 1..n_cat-1
  # cum_probs is n x (n_cat - 1) matrix
  cum_probs <- matrix(nrow = n, ncol = n_cat - 1L)
  for (k in seq_len(n_cat - 1L)) {
    cum_probs[, k] <- stats::plogis(eta + d_vec[k])
  }
  # Response (0-indexed) = number of thresholds exceeded
  u <- stats::runif(n)
  as.integer(rowSums(cum_probs > u))
}


#' Generate theta values from age-informed prior
#'
#' Regression coefficients estimated from the NE25 validation sample
#' (N = 2,785) using the scoring model's design matrix
#' \code{[1, ln(age + 1), age]}.
#'
#' @param years_old Numeric vector of ages.
#' @param n Integer number of observations.
#' @param factor_names Character vector of factor names.
#' @return Numeric matrix with n rows and 7 columns.
#' @noRd
generate_theta_from_prior <- function(years_old, n, factor_names) {
  # NE25 regression: rows are F/GEN/EAT/EXT/INT/SLE/SOC
  # columns are intercept, ln(age+1), age
  beta <- matrix(
    c(
      -8.7743, 9.7330, -1.0462,
      -1.3284, 2.0275, -0.4251,
      -0.2053, 0.7138, -0.1733,
      -0.0218, 0.6678, -0.1855,
      -0.1615, 0.8208, -0.2070,
      0.1845, 0.3143, -0.1450,
      0.0060, 0.6241, -0.1527
    ),
    nrow = 7, ncol = 3, byrow = TRUE
  )
  sigma <- c(1.194, 1.073, 0.757, 0.903, 1.050, 0.778, 0.765)

  # Design matrix: [1, ln(age + 1), age]
  x <- cbind(1, log(years_old + 1), years_old)

  # mu = X %*% t(beta) -> n x 7 matrix
  mu <- x %*% t(beta)

  # Draw theta ~ N(mu, sigma)
  theta <- matrix(nrow = n, ncol = 7)
  for (k in seq_len(7)) {
    theta[, k] <- stats::rnorm(n, mean = mu[, k], sd = sigma[k])
  }
  colnames(theta) <- factor_names
  theta
}
