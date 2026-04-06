// grm_unidimensional.stan
// Graded Response Model — single factor (F), MAP scoring via BFGS
//
// Fixed item parameters (slopes + thresholds from Mplus calibration).
// Age-informed prior on theta via QR-reparameterized latent regression.
// FIML: missing responses skipped via precomputed observed-item indices.
//
// Thresholds are sign-flipped from Mplus convention: d = -tau.
// GRM probability: P(Y >= k | theta) = inv_logit(a * theta + d[k])

functions {
  /**
   * Log-probability of a GRM response.
   *
   * @param y      Observed response (1-indexed: 1, ..., n_cat)
   * @param eta    Linear predictor: a * theta
   * @param d_row  Row of sign-flipped thresholds (padded; uses indices 1..n_cat-1)
   * @param n_cat  Number of response categories for this item
   * @return       Log-probability of observing y
   */
  real grm_log_prob(int y, real eta, row_vector d_row, int n_cat) {
    if (y == 1) {
      return log1m_inv_logit(eta + d_row[1]);
    } else if (y == n_cat) {
      return log_inv_logit(eta + d_row[n_cat - 1]);
    } else {
      return log_diff_exp(
        log_inv_logit(eta + d_row[y - 1]),
        log_inv_logit(eta + d_row[y])
      );
    }
  }
}

data {
  int<lower=1> N;                       // number of persons
  int<lower=1> J;                       // number of items
  int<lower=2> max_cat;                 // max response categories across items
  array[J] int<lower=2> n_cat;         // number of categories per item
  array[N, J] int<lower=0> y;          // responses: 0 = missing, 1..n_cat[j] = observed
  vector[J] a;                          // fixed item slopes
  matrix[J, max_cat - 1] d;            // fixed thresholds (d = -tau, padded with 0)
  matrix[N, 3] X;                       // design matrix: [1, ln(age + 1), age]
}

transformed data {
  // QR decomposition of design matrix for numerical stability
  matrix[N, 3] Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  matrix[3, 3] R_ast = qr_thin_R(X) / sqrt(N - 1);
  matrix[3, 3] R_ast_inverse = inverse(R_ast);

  // FIML: precompute observed item indices per person
  array[N] int<lower=0> n_obs;
  array[N, J] int obs_idx;
  for (i in 1:N) {
    n_obs[i] = 0;
    for (j in 1:J) {
      obs_idx[i, j] = 0;
      if (y[i, j] > 0) {
        n_obs[i] += 1;
        obs_idx[i, n_obs[i]] = j;
      }
    }
  }
}

parameters {
  vector[N] theta;                      // latent trait (MAP estimates)
  vector[3] gamma;                      // latent regression coefficients (QR space)
}

model {
  // Prior on regression coefficients
  gamma ~ normal(0, 2.5);

  // Age-informed prior: theta_i ~ N(X_i * beta, 1)
  theta ~ normal(Q_ast * gamma, 1);

  // FIML likelihood: loop over persons and their observed items
  for (i in 1:N) {
    for (k in 1:n_obs[i]) {
      int j = obs_idx[i, k];
      real eta = a[j] * theta[i];
      target += grm_log_prob(y[i, j], eta, d[j], n_cat[j]);
    }
  }
}

generated quantities {
  // Back-transform regression coefficients to original scale
  vector[3] beta = R_ast_inverse * gamma;
}
