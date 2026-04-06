test_that("score_kidsights errors on missing id_cols", {
  df <- data.frame(AA4 = c(0, 1), years_old = c(2, 3))
  expect_error(score_kidsights(df, id_cols = "pid"), "ID columns not found")
})

test_that("score_kidsights errors on missing years_old", {
  df <- data.frame(AA4 = c(0, 1), pid = c(1, 2))
  expect_error(score_kidsights(df, id_cols = "pid"), "years_old")
})

test_that("score_kidsights errors on duplicate IDs", {
  df <- data.frame(AA4 = c(0, 1), years_old = c(2, 3), pid = c(1, 1))
  expect_error(score_kidsights(df, id_cols = "pid"), "uniquely identify")
})

test_that("score_kidsights errors on negative age", {
  df <- data.frame(AA4 = c(0, 1), years_old = c(-1, 3), pid = c(1, 2))
  expect_error(score_kidsights(df, id_cols = "pid"), "negative")
})

# Helper to find validation data (works both interactively and during R CMD check)
find_validation_file <- function(filename) {
  paths <- c(
    test_path("..", "..", "data-raw", "validation", filename),
    file.path("data-raw", "validation", filename)
  )
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  NULL
}

test_that("score_kidsights supports multiple id_cols", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  items_path <- find_validation_file("validation_ne25_items.rda")
  ages_path <- find_validation_file("validation_ne25_ages.rda")
  skip_if(is.null(items_path), "Validation data not available")

  load(items_path)
  load(ages_path)

  df <- validation_ne25_items[1:50, ]
  df$years_old <- validation_ne25_ages$years_old[1:50]
  df$site <- "ne25"
  df$rid <- seq_len(50)

  result <- score_kidsights(df, id_cols = c("site", "rid"))
  expect_true("site" %in% names(result))
  expect_true("rid" %in% names(result))
  expect_true("theta" %in% names(result))
  expect_equal(nrow(result), 50)
})

test_that("score_kidsights returns NA for observations below min_responses", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  items_path <- find_validation_file("validation_ne25_items.rda")
  ages_path <- find_validation_file("validation_ne25_ages.rda")
  skip_if(is.null(items_path), "Validation data not available")

  load(items_path)
  load(ages_path)

  df <- validation_ne25_items[1:20, ]
  df$years_old <- validation_ne25_ages$years_old[1:20]
  df$rid <- seq_len(20)

  item_cols <- setdiff(names(df), c("years_old", "rid"))
  df[1, item_cols] <- NA

  result <- score_kidsights(df, id_cols = "rid", min_responses = 5)
  expect_true(is.na(result$theta[1]))
  expect_equal(nrow(result), 20)
})

test_that("score_kidsights warns on out-of-range responses", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  items_path <- find_validation_file("validation_ne25_items.rda")
  ages_path <- find_validation_file("validation_ne25_ages.rda")
  skip_if(is.null(items_path), "Validation data not available")

  load(items_path)
  load(ages_path)

  df <- validation_ne25_items[1:50, ]
  df$years_old <- validation_ne25_ages$years_old[1:50]
  df$rid <- seq_len(50)

  df$AA4[1] <- 99

  expect_warning(
    score_kidsights(df, id_cols = "rid"),
    "outside valid range"
  )
})
