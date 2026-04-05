# KidsightsPublic Package Rebuild: CmdStan MAP Scoring

**Created:** 2026-04-04 | **Status:** Draft | **Repo:** C:/Users/marcu/git-repositories/KidsightsPublic

---

## Context

The KidsightsPublic R package currently uses `mirt` for IRT person scoring with fixed item parameters from Mplus calibration. We're rebuilding it to use **CmdStan BFGS** for MAP scoring, eliminating the Mplus dependency and enabling automated scoring for any study site.

**Current state:** 2 functions (`kidsights_scores()`, `early_learning()`), mirt-based, no tests.
**Target:** CmdStan-based MAP scorer with auto-lexicon detection, GRM support, age-informed priors, and validated against Mplus output.

---

## Model Specification

### Item Response Model: Graded Response Model (GRM)

Items have varying numbers of response categories (2 to 6). For item i with C_i categories:

```
P(Y_i >= c | theta) = logistic(a_i * lambda_{i,k} * theta_k + d_{i,c})
```

where `a_i` = slope, `lambda_{i,k}` = factor loading pattern (0/1), `d_{i,c}` = threshold c.

All item parameters (`a`, `d`) are **fixed data** — not estimated.

### Latent Regression Prior (Age-Informed)

Each factor k has its own age-informed prior:

```
theta_{i,k} ~ N(X_i * beta_k, 1)
```

where `X_i = [1, ln(age_i + 1), age_i]` (intercept + 2 basis functions).

- **QR decomposition** on the N × 3 design matrix X for numerical stability
- **Residual variance fixed at 1** (standard normal)
- **Betas estimated per dataset** (jointly optimized with all thetas)

### Two Models

1. **Unidimensional** (Kidsights overall): 1 factor, 220 items, N+3 parameters
2. **Bifactor** (psychosocial): 7 factors (general + eat/ext/int/sle/soc/gen_ps), 220 items, 7N+21 parameters

### Optimization

- **Method:** BFGS via CmdStan `$optimize()`
- **Output:** MAP theta estimates + Hessian → CSEMs via inverse Fisher information

---

## Data Wrangling Plan

### Phase 1: Extract Item Parameters from Mplus .out Files

**Source files:**
- `Kidsights-Data-Platform/calibration/ne25/manual_2023_scale/mplus/kidsights_2023_calibration_ne25.out` (unidimensional)
- `Kidsights-Data-Platform/calibration/ne25/manual_2023_scale/mplus/all_2023_calibration_ne25.out` (bifactor)

**Parse and extract:**
- Item names (220 items)
- Slopes (a) per item per factor — fixed (@) and freely estimated values
- Thresholds (d) per item — $1 through $5 (varying by item)
- Factor loading pattern (which items load on which factors in bifactor model)
- Number of categories per item (determined by max threshold index)

**Output format:** R list or data frame stored as package data:
```r
item_params$unidimensional  # slopes + thresholds for 220 items × 1 factor
item_params$bifactor        # slopes + thresholds for 220 items × 7 factors
```

**Item category counts (from Mplus .out):**
- ~220 items with $1 (all items have at least 1 threshold)
- ~44 items with $2 (3+ categories)
- ~42 items with $3 (4+ categories)
- ~33 items with $4 (5+ categories)
- ~1 item with $5 (6 categories)

### Phase 2: Build Codebook Lexicon Lookup

**Source:** `Kidsights-Data-Platform/codebook/data/codebook.json`

Bundle a lookup table that maps column names → canonical item IDs for any lexicon:
- `C020` (ne25/mn26) → `AA4`
- `crosec004` (gsed) → `AA4`
- `LF5` (credi) → `AA4`
- `AA4` (equate) → `AA4`

**Auto-detection:** Given a data frame, check column names against all lexicons, find the best match (highest overlap), and return the lexicon name + column mapping.

### Phase 3: Write Stan Models

**Key performance decisions:**
1. **Missing data via `transformed data` indexing** — precompute observed item indices per person once, avoid branching in likelihood
2. **Vectorized `ordered_logistic_lpmf`** — push inner item loop into Stan's C++ implementation
3. **`reduce_sum` parallelism** — parallelize person-level likelihood across CPU threads
4. **QR decomposition** — reparameterize design matrix for numerical stability
5. **Threshold sign flip** — Mplus threshold τ converted to Stan intercept d = -τ

**CRITICAL: Mplus parameterizes as `a*theta - tau`. Stan's `ordered_logistic` uses `a*theta + d`. Therefore `d = -tau_mplus` (sign flip on all thresholds when parsing .out files).**

**`inst/stan/grm_unidimensional.stan` (sketch):**
```stan
functions {
  // Per-person log-likelihood (for reduce_sum parallelism)
  real person_lpdf(array[] int y_slice,
                   int start, int end,
                   vector theta, vector a, matrix d,
                   array[,] int obs_idx, array[] int n_obs,
                   array[] int n_cat) {
    real ll = 0;
    for (i in start:end) {
      // Vectorized: extract observed items for person i
      int n = n_obs[i];
      // Use ordered_logistic_lpmf on observed items only
      for (k in 1:n) {
        int j = obs_idx[i, k];
        // ordered_logistic with cutpoints = -d (Stan convention)
        ll += ordered_logistic_lpmf(y_slice[i - start + 1, k] |
                                     a[j] * theta[i],
                                     segment(d[j], 1, n_cat[j] - 1));
      }
    }
    return ll;
  }
}

data {
  int<lower=1> N;                    // persons
  int<lower=1> J;                    // items
  int<lower=1> max_cat;              // max categories across items
  array[J] int<lower=2> n_cat;      // categories per item
  array[N, J] int y;                // responses (-1 = missing)
  vector[J] a;                       // fixed slopes
  matrix[J, max_cat-1] d;           // fixed thresholds (sign-flipped from Mplus)
  matrix[N, 3] X;                    // design matrix [1, ln(age+1), age]
}

transformed data {
  // QR decomposition of design matrix
  matrix[N, 3] Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  matrix[3, 3] R_ast = qr_thin_R(X) / sqrt(N - 1);
  matrix[3, 3] R_ast_inverse = inverse(R_ast);

  // Precompute observed item indices per person (FIML)
  array[N] int n_obs;
  array[N, J] int obs_idx;
  for (i in 1:N) {
    n_obs[i] = 0;
    for (j in 1:J) {
      if (y[i, j] >= 0) {
        n_obs[i] += 1;
        obs_idx[i, n_obs[i]] = j;
      }
    }
  }
}

parameters {
  vector[N] theta;                   // person MAP scores
  vector[3] gamma;                   // latent regression (QR space)
}

model {
  // Prior on regression coefficients
  gamma ~ normal(0, 2.5);

  // Age-informed prior on theta: theta_i ~ N(X_i * beta, 1)
  theta ~ normal(Q_ast * gamma, 1);

  // FIML likelihood: only observed items per person
  // Use reduce_sum for parallelism across persons
  target += reduce_sum(person_lpdf, y, 1,
                       theta, a, d, obs_idx, n_obs, n_cat);
}

generated quantities {
  // Back-transform regression coefficients
  vector[3] beta = R_ast_inverse * gamma;
}
```

**`inst/stan/grm_bifactor.stan`:**
Same structure extended to K=7 factors:
- `matrix[N, K] theta` (7 scores per person)
- `matrix[3, K] gamma` (regression coefficients per factor, QR space)
- `matrix[J, K] lambda` (0/1 factor loading pattern: which items load on which factors)
- Item likelihood: `sum_k(a[j] * lambda[j,k] * theta[i,k]) + d[j,c]`
- Each factor gets own age-informed prior: `theta[i,k] ~ N(Q_ast[i] * gamma[,k], 1)`
- `reduce_sum` parallelizes across persons (each person's 7D theta is independent given beta)

### Phase 4: Build R Functions

**Core functions:**
- `parse_mplus_out(filepath)` → Extract item params from .out file
- `detect_lexicon(data, codebook)` → Auto-detect which lexicon column names match
- `map_to_canonical(data, lexicon)` → Rename columns to canonical equate IDs
- `prepare_stan_data(data, item_params, model_type)` → Build Stan data list
- `score_kidsights(data, model = "unidimensional", ...)` → Main user-facing function
- `score_psychosocial(data, ...)` → Bifactor scoring

**User-facing API:**
```r
# Simplest usage — auto-detects lexicon, scores all applicable domains
scores <- score_kidsights(data)

# With options
scores <- score_kidsights(data, model = "bifactor", min_responses = 5)
```

### Phase 5: Test Data Pipeline + Validation

#### Test Data Preparation

**Source:** NE25 pipeline data from Kidsights-Data-Platform repo.

**Script:** `data-raw/prepare_test_data.R` — pulls test data from the Platform and formats for package testing:

1. Load `mplus_dat.dat` with column names from `mplus_dat.inp` (2,785 records, 233 columns)
   - Source: `Kidsights-Data-Platform/calibration/ne25/manual_2023_scale/mplus/`
   - Columns: rid, pid, recordid, covariates (logyrs, yrs3, ...), 220 item responses (AA4-PS049)
   - Items already in equate lexicon (AA4, BB5, CC4, etc.)
   - Missing coded as `.` (Mplus convention) → convert to NA

2. Extract `years_old` from `ne25_transformed` table in DuckDB (or back-compute from `logyrs`)
   - **Age must be in raw years** — Stan model computes `ln(years_old + 1)` and `years_old` internally
   - Design matrix: `X = [1, ln(years_old + 1), years_old]` (3 columns)

3. Load Mplus ground truth scores:
   - `ne25_scores_kidsights_2023_scale.dat` (2,781 records, 237 columns — last 2 are theta + SE)
   - `ne25_scores_all_2023_scale.dat` (2,785 records, 293 columns — last 14 are 7 thetas + 7 SEs)

4. Save as `data/test_ne25_items.rda`, `data/test_ne25_ages.rda`, `data/test_mplus_scores_unidim.rda`, `data/test_mplus_scores_bifactor.rda`

#### Validation Strategy

**Important:** The Mplus model used **11 covariates** in the latent regression:
```
F ON logyrs yrs3 school logfpl phq2 schXyrs3 fplXyrs3 phqXyrs3 black hisp other;
```

The new CmdStan model uses only **3 basis functions** of age: `[1, ln(years_old + 1), years_old]`.

This means scores will NOT be identical to Mplus — the priors differ. Two-stage validation:

**Validation: Approximate match against Mplus**
- Run CmdStan with production prior `[1, ln(years_old + 1), years_old]`
- Compare theta estimates against Mplus scores (which used 11 covariates including FPL, PHQ2, race)
- Expected: correlation > 0.98 (likelihood dominates prior for most children)
- Scores will diverge more for children with few item responses (prior matters more)
- Verify: age-score gradient is positive, scores are in reasonable range, CSEMs are plausible
- This is NOT an exact replication — different priors, same GRM likelihood + item parameters

---

## Implementation Sprints

### Sprint 1: Mplus Parser + Item Parameters + Test Data
- `R/parse_mplus_out.R` — parser for .out files (slopes, thresholds with sign flip d = -tau)
- `data-raw/extract_params.R` — script to create item_params from .out files
- `data-raw/prepare_test_data.R` — pull NE25 test data from Kidsights-Data-Platform
- `data/item_params.rda` — bundled item parameters (220 items, unidim + bifactor)
- `data/codebook_lookup.rda` — lexicon lookup table from codebook.json
- `data/test_ne25_items.rda` — NE25 item responses (equate lexicon, 2785 records)
- `data/test_ne25_ages.rda` — NE25 ages in years (raw, not log-transformed)
- `data/test_mplus_scores_unidim.rda` — Mplus theta + SE ground truth
- `data/test_mplus_scores_bifactor.rda` — Mplus 7-factor scores ground truth
- Test: Extracted params match .out file values; test data loads correctly

### Sprint 2: Stan Models
- `inst/stan/grm_unidimensional.stan` — single factor GRM with age prior + QR
- `inst/stan/grm_bifactor.stan` — 7-factor bifactor GRM with age priors + QR
- Test: Models compile, optimize on synthetic data

### Sprint 3: R Scoring Functions
- `R/detect_lexicon.R` — auto-detection
- `R/score_kidsights.R` — main scoring function (replaces mirt version)
- `R/score_psychosocial.R` — bifactor scoring
- Update DESCRIPTION: remove mirt, add cmdstanr
- Test: Runs on NE25 data, produces scores

### Sprint 4: Validation + Testing
- `tests/testthat/test-score_kidsights.R` — equivalence to Mplus
- `tests/testthat/test-detect_lexicon.R` — lexicon detection
- `tests/testthat/test-parse_mplus.R` — parser correctness
- Validation report: correlation, RMSE, scatter plots
- Documentation: vignette with usage examples

### Sprint 5: Package Polish
- Update DESCRIPTION, README, NAMESPACE
- Add `early_learning()` replacement (or deprecate if subsumed by bifactor)
- roxygen2 documentation for all exported functions
- `pkgdown` site or vignettes

---

## Files to Create

| File | Sprint | Purpose |
|------|--------|---------|
| `todo/scoring_package_plan.md` | 0 | This plan |
| `R/parse_mplus_out.R` | 1 | Mplus .out parser |
| `data-raw/extract_params.R` | 1 | Parameter extraction script |
| `data-raw/prepare_test_data.R` | 1 | Pull NE25 test data from Platform repo |
| `inst/stan/grm_unidimensional.stan` | 2 | Unidimensional GRM |
| `inst/stan/grm_bifactor.stan` | 2 | Bifactor GRM |
| `R/detect_lexicon.R` | 3 | Auto-lexicon detection |
| `R/score_kidsights.R` | 3 | Main scoring (replaces mirt version) |
| `R/score_psychosocial.R` | 3 | Bifactor scoring |
| `tests/testthat/test-*.R` | 4 | Test suite |

## Files to Modify

| File | Sprint | Change |
|------|--------|--------|
| `DESCRIPTION` | 3 | Remove mirt, add cmdstanr, update metadata |
| `NAMESPACE` | 3 | Update exports |
| `README.md` | 5 | Rewrite for new API |

## Files to Remove (or Archive)

| File | Reason |
|------|--------|
| `R/kidsights_scores.R` | Replaced by CmdStan version |
| `R/early_learning.R` | Replaced or deprecated |
| `data/internals.rda` | Replaced by parsed Mplus params |
| `data/calibdat.rda` | May keep for validation, remove from package |

---

## Verification

- **Sprint 1:** Extracted params match manual spot-check of .out files; test data loads; ages in years (not log)
- **Sprint 2:** Stan models compile and optimize on 100-person synthetic data
- **Sprint 3:** Scores produced for NE25 test data, reasonable range
- **Sprint 4:** Production model `[1, ln(age+1), age]` — correlation > 0.98 with Mplus, age gradient positive
- **Sprint 5:** `R CMD check` passes, package installs cleanly

**Note on age:** All age values are in **raw years** (e.g., 2.5 = 2 years 6 months). The Stan model computes `ln(years_old + 1)` internally. R prep functions should NOT pre-transform age.

## Risk Register

| Risk | Mitigation |
|------|-----------|
| BFGS convergence issues for large N × 7K | Start with unidimensional, tune initial values from marginal estimates |
| GRM likelihood gradient complexity | Stan handles autodiff; verify against finite differences |
| Mplus .out parsing fragile | Multiple test files, structured regex patterns |
| CmdStan installation dependency | Document setup, consider fallback to rstan |
| MAP ≠ EAP (Mplus uses different estimator) | Set tolerance empirically, document expected differences |
