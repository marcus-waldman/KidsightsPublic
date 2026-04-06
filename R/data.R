#' Item parameters from Mplus GRM calibration
#'
#' Fixed item parameters extracted from the Mplus bifactor calibration
#' (\code{all_2023_calibration_ne25.out}) using \code{MplusAutomation::readModels()}.
#' Contains slopes, thresholds, category counts, and factor names for 265 items
#' across 7 factors.
#'
#' @format A list with four components:
#' \describe{
#'   \item{slopes}{Data frame with columns \code{item}, \code{factor}, and
#'     \code{slope}. 314 rows covering 7 factors (F, GEN, EAT, EXT, INT, SLE,
#'     SOC). All slopes are fixed (from calibration).}
#'   \item{thresholds}{Data frame with columns \code{item}, \code{threshold}
#'     (integer index), and \code{d} (sign-flipped: \code{d = -tau_mplus}).
#'     430 rows.}
#'   \item{n_cat}{Named integer vector of response category counts per item
#'     (range 2-6). Length 265.}
#'   \item{factors}{Character vector of factor names:
#'     \code{c("EAT", "EXT", "F", "GEN", "INT", "SLE", "SOC")}.}
#' }
#'
#' @details
#' The 220 non-PS items load on factor F (Kidsights overall). The 44 PS items
#' load on GEN (general psychosocial) plus one subfactor (EAT, EXT, INT, SLE,
#' or SOC) but not on F. One item (CC79Y) has a case variant between Mplus
#' (all-caps) and the codebook (mixed case).
#'
#' Thresholds are stored as \code{d = -tau} where \code{tau} is the Mplus
#' parameterization. The GRM probability is
#' \code{P(Y >= k | theta) = inv_logit(a * theta + d[k])}.
#'
#' @source Mplus calibration: \code{all_2023_calibration_ne25.out}
"item_params"


#' Codebook lexicon lookup table
#'
#' Maps canonical equate item IDs to site-specific column names across
#' multiple lexicons. Used by \code{\link{detect_lexicon}} to auto-identify
#' which lexicon a dataset uses.
#'
#' @format A data frame with 1802 rows and 3 columns:
#' \describe{
#'   \item{equate_id}{Canonical item identifier (e.g., \code{"AA4"}).}
#'   \item{lexicon}{Lexicon name (e.g., \code{"ne25"}, \code{"credi"},
#'     \code{"gsed"}).}
#'   \item{column_name}{Column name used in that lexicon (e.g., \code{"C020"}
#'     for AA4 in the ne25 lexicon).}
#' }
#'
#' @details
#' Available lexicons: cahmi21, cahmi22, credi, ecdi, equate, gsed, kidsight,
#' mn26, ne20, ne22, ne25. The equate lexicon uses the canonical IDs as
#' column names.
#'
#' @source \code{Kidsights-Data-Platform/codebook/data/codebook.json}
"codebook_lookup"
