// grm_bifactor.stan
// Graded Response Model — bifactor (GEN + 5 subfactors), MAP scoring via BFGS
//
// 44 PS items, each loading on GEN plus exactly one subfactor (EAT/EXT/INT/SLE/SOC).
// Fixed item parameters (slopes + thresholds from Mplus calibration).
// Per-factor age-informed priors via QR-reparameterized latent regression.
// FIML: missing responses skipped via precomputed observed-item indices.
//
// Linear predictor: eta = a_gen[j] * theta_gen[i] + a_sub[j] * theta_sub[i, sub_idx[j]]

functions {
  /**
   * Log-probability of a GRM response.
   *
   * @param y      Observed response (1-indexed: 1, ..., n_cat)
   * @param eta    Linear predictor
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
  int<lower=1> J;                       // number of items (44 PS items)
  int<lower=1> K_sub;                   // number of subfactors (5: EAT, EXT, INT, SLE, SOC)
  int<lower=2> max_cat;                 // max response categories across items
  array[J] int<lower=2> n_cat;         // number of categories per item
  array[N, J] int<lower=0> y;          // responses: 0 = missing, 1..n_cat[j] = observed
  vector[J] a_gen;                      // fixed GEN slopes per item
  vector[J] a_sub;                      // fixed subfactor slopes per item
  array[J] int<lower=1, upper=K_sub> sub_idx;  // which subfactor each item loads on
  matrix[J, max_cat - 1] d;            // fixed thresholds (d = -tau, padded with 0)
  matrix[N, 3] X;                       // design matrix: [1, ln(age + 1), age]
}

transformed data {
  int K = K_sub + 1;                    // total factors (GEN + subfactors)

  // QR decomposition of design matrix
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
  vector[N] theta_gen;                  // GEN factor scores
  matrix[N, K_sub] theta_sub;          // subfactor scores (EAT, EXT, INT, SLE, SOC)
  matrix[3, K] gamma;                   // regression coefficients (QR space), col order: GEN, sub1..sub5
}

model {
  // Priors on regression coefficients
  to_vector(gamma) ~ normal(0, 2.5);

  // Per-factor age-informed priors
  // GEN: theta_gen ~ N(Q_ast * gamma[, 1], 1)
  theta_gen ~ normal(Q_ast * gamma[, 1], 1);

  // Subfactors: theta_sub[, k] ~ N(Q_ast * gamma[, k+1], 1)
  for (k in 1:K_sub) {
    theta_sub[, k] ~ normal(Q_ast * gamma[, k + 1], 1);
  }

  // FIML likelihood
  for (i in 1:N) {
    for (m in 1:n_obs[i]) {
      int j = obs_idx[i, m];
      real eta = a_gen[j] * theta_gen[i]
               + a_sub[j] * theta_sub[i, sub_idx[j]];
      target += grm_log_prob(y[i, j], eta, d[j], n_cat[j]);
    }
  }
}

generated quantities {
  // Back-transform regression coefficients to original scale
  matrix[3, K] beta;
  for (k in 1:K) {
    beta[, k] = R_ast_inverse * gamma[, k];
  }
}
