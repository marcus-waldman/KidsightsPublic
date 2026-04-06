# KidsightsPublic

R package for computing Kidsights developmental screening scores using Graded
Response Models (GRM) with fixed item parameters. Scores are estimated via MAP
(maximum a posteriori) optimization using CmdStan.

## Installation

### Prerequisites

CmdStan must be installed first:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

### Install the package

```r
devtools::install_github("marcus-waldman/KidsightsPublic")
```

## Usage

Your data frame needs:

- Item response columns (0-indexed) with column names from a recognized lexicon
- A `years_old` column with age in years
- One or more ID columns

### Kidsights overall score

```r
library(KidsightsPublic)

results <- score_kidsights(data, id_cols = c("child_id", "site"))
head(results)
#>   child_id site  theta
#>   1001     ne25  -2.31
#>   1002     ne25   3.45
```

### Psychosocial domain scores

```r
ps_results <- score_psychosocial(data, id_cols = "child_id")
head(ps_results)
#>   child_id theta_gen theta_eat theta_ext theta_int theta_sle theta_soc
#>   1001      0.12     -0.45      0.33     -0.21      0.55     -0.10
```

### Lexicon auto-detection

The package auto-detects which lexicon your column names follow (equate,
ne25, credi, gsed, mn26, and others):

```r
detect_lexicon(data)
#> $lexicon
#> [1] "ne25"
#> $n_matched
#> [1] 185
```

## Scoring details

- **Model**: GRM with fixed item parameters from Mplus calibration (NE25)
- **Optimization**: L-BFGS via CmdStan (typically 2-15 seconds)
- **Prior**: Age-informed latent regression: `theta ~ N(beta_0 + beta_1 * ln(age + 1) + beta_2 * age, 1)`
- **Missing data**: FIML (items with missing responses are skipped in the likelihood)
- **Minimum responses**: Observations with fewer than 5 item responses are not scored (configurable)

## License

MIT
