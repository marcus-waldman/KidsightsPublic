test_that("detect_lexicon identifies equate lexicon", {
  df <- data.frame(AA4 = c(0, 1), AA5 = c(1, 0), AA6 = c(0, 1))
  result <- detect_lexicon(df, min_overlap = 0)
  expect_equal(result$lexicon, "equate")
  expect_true(result$n_matched >= 3)
})

test_that("detect_lexicon is case-insensitive", {
  df <- data.frame(aa4 = c(0, 1), aa5 = c(1, 0), Aa6 = c(0, 1))
  result <- detect_lexicon(df, min_overlap = 0)
  expect_equal(result$lexicon, "equate")
})

test_that("detect_lexicon identifies ne25 lexicon", {
  df <- data.frame(
    C020 = c(0, 1), C023 = c(1, 0), C026 = c(0, 1),
    C033 = c(1, 1), C038 = c(0, 0)
  )
  result <- detect_lexicon(df, min_overlap = 0)
  expect_equal(result$lexicon, "ne25")
  expect_true("AA4" %in% result$mapping)
})

test_that("detect_lexicon returns correct mapping structure", {
  df <- data.frame(AA4 = c(0, 1), AA5 = c(1, 0))
  result <- detect_lexicon(df, min_overlap = 0)
  expect_type(result$mapping, "character")
  expect_true(all(names(result$mapping) %in% names(df)))
  expect_true("AA4" %in% result$mapping)
  expect_true("AA5" %in% result$mapping)
})

test_that("detect_lexicon errors on unrecognized columns", {
  df <- data.frame(xyz1 = c(0, 1), xyz2 = c(1, 0))
  expect_error(detect_lexicon(df), "No lexicon matched")
})
