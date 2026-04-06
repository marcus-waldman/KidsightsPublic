# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

**KidsightsPublic** is an R package for IRT-based person scoring using the Graded Response Model (GRM) with fixed item parameters calibrated in Mplus. Currently being rebuilt from `mirt`-based scoring to **CmdStan BFGS MAP scoring**.

- **Primary Users**: Research teams at study sites who need automated developmental screening scores
- **Plan**: See `todo/scoring_package_plan.md` for full implementation plan

## Development Workflow

Before committing any code, run this sequence:

1. `devtools::load_all()` — load package for testing
2. `devtools::document()` — update documentation from roxygen2
3. `lintr::lint_package()` — check Tidyverse style compliance
4. `styler::style_pkg()` — auto-fix formatting
5. `devtools::test()` — run all tests
6. `devtools::check()` — full R CMD check (must pass with 0 errors, 0 warnings, 0 notes)

## Style Requirements

**Strictly follows the [Tidyverse Style Guide](https://style.tidyverse.org/)**

- **snake_case** naming for functions and variables
- **120 character** line length maximum
- **2-space indentation** (never tabs)
- **roxygen2 documentation** required for all exported functions
- **No `browser()`, `print()`, or `library()` calls** in R/ source files
- Use `@importFrom pkg fun` or `pkg::fun()` for external functions

## CRAN Compliance

The package must always be CRAN-submission ready:

- All exported functions have complete roxygen2 documentation with `@param`, `@return`, `@export`, `@examples`
- No absolute file paths in package code (use `system.file()` for package files)
- All dependencies declared in DESCRIPTION Imports/Suggests
- Examples must run (or use `\donttest{}` / `\dontrun{}`)
- No `library()` or `require()` calls inside functions — use `::` or `@importFrom`

## Key Technical Notes

- **Threshold sign flip**: Mplus parameterizes as `a*theta - tau`; Stan uses `a*theta + d`. Always apply `d = -tau` when converting from Mplus.
- **Age in raw years**: All age values stored as raw years (e.g., 2.5). Stan model computes `ln(years_old + 1)` internally. R functions should NOT pre-transform age.
- **MplusAutomation**: Use `MplusAutomation::readModels()` for parsing Mplus .out files — do not write custom parsers.
