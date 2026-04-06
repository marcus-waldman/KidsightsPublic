test_that("item_params has correct structure", {
  expect_type(item_params, "list")
  expect_named(item_params, c("slopes", "thresholds", "n_cat", "factors"))
})

test_that("item_params has correct counts", {
  expect_equal(nrow(item_params$slopes), 314)
  expect_equal(nrow(item_params$thresholds), 430)
  expect_equal(length(item_params$n_cat), 265)
  expect_equal(length(item_params$factors), 7)
})

test_that("item_params factors are correct", {
  expect_equal(
    sort(item_params$factors),
    c("EAT", "EXT", "F", "GEN", "INT", "SLE", "SOC")
  )
})

test_that("item_params slopes have correct columns", {
  expect_named(item_params$slopes, c("item", "factor", "slope"))
})

test_that("item_params thresholds are sign-flipped from Mplus", {
  # AA4$1 in Mplus: tau = -5.920, so d = -(-5.920) = 5.920
  aa4_d1 <- item_params$thresholds$d[
    item_params$thresholds$item == "AA4" & item_params$thresholds$threshold == 1
  ]
  expect_equal(aa4_d1, 5.920)
})

test_that("item_params slopes match known values", {
  # AA4 F slope = 0.430
  aa4_slope <- item_params$slopes$slope[
    item_params$slopes$item == "AA4" & item_params$slopes$factor == "F"
  ]
  expect_equal(aa4_slope, 0.430)

  # PS003 EXT slope = 2.069
  ps003_slope <- item_params$slopes$slope[
    item_params$slopes$item == "PS003" & item_params$slopes$factor == "EXT"
  ]
  expect_equal(ps003_slope, 2.069)
})

test_that("item_params n_cat values are plausible", {
  expect_true(all(item_params$n_cat >= 2))
  expect_true(all(item_params$n_cat <= 6))
})

test_that("F factor has 220 items", {
  f_count <- sum(item_params$slopes$factor == "F")
  expect_equal(f_count, 220)
})

test_that("GEN factor has 44 items", {
  gen_count <- sum(item_params$slopes$factor == "GEN")
  expect_equal(gen_count, 44)
})
