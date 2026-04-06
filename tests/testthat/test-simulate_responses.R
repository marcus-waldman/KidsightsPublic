test_that("simulate_responses returns correct structure", {
  sim <- simulate_responses(years_old = 3, n = 5, seed = 42)

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 5)
  # id + years_old + 264 items = 266 columns
  expect_equal(ncol(sim), 266)
  expect_true("id" %in% names(sim))
  expect_true("years_old" %in% names(sim))
  expect_equal(sim$id, 1:5)
  expect_equal(sim$years_old, rep(3, 5))
})

test_that("simulate_responses has correct item count", {
  sim <- simulate_responses(years_old = 3, seed = 1)

  ip <- KidsightsPublic::item_params
  f_items <- ip$slopes$item[ip$slopes$factor == "F"]
  subfactor_items <- unique(ip$slopes$item[ip$slopes$factor != "F"])
  f_only <- setdiff(f_items, subfactor_items)
  ps_items <- ip$slopes$item[ip$slopes$factor == "GEN"]

  item_cols <- setdiff(names(sim), c("id", "years_old"))
  expect_equal(length(item_cols), length(f_only) + length(ps_items))
  expect_true(all(f_only %in% item_cols))
  expect_true(all(ps_items %in% item_cols))
})

test_that("simulate_responses produces valid response ranges", {
  sim <- simulate_responses(years_old = 3, n = 50, seed = 99)
  ip <- KidsightsPublic::item_params

  item_cols <- setdiff(names(sim), c("id", "years_old"))
  for (item in item_cols) {
    vals <- sim[[item]]
    expect_true(all(!is.na(vals)), info = paste("NA found in", item))
    expect_true(all(vals >= 0), info = paste("Negative value in", item))
    expect_true(
      all(vals <= ip$n_cat[item] - 1L),
      info = paste("Value exceeds max category in", item)
    )
  }
})

test_that("simulate_responses is reproducible with seed", {
  sim1 <- simulate_responses(years_old = 2, n = 10, seed = 123)
  sim2 <- simulate_responses(years_old = 2, n = 10, seed = 123)
  expect_identical(sim1, sim2)
})

test_that("simulate_responses seed does not affect caller RNG", {
  set.seed(999)
  before <- stats::runif(1)
  set.seed(999)
  simulate_responses(years_old = 3, n = 5, seed = 42)
  after <- stats::runif(1)
  expect_equal(before, after)
})

test_that("simulate_responses handles vector ages", {
  ages <- c(1, 2.5, 4, 5.5)
  sim <- simulate_responses(years_old = ages, seed = 10)
  expect_equal(nrow(sim), 4)
  expect_equal(sim$years_old, ages)
})

test_that("simulate_responses accepts custom theta", {
  theta_mat <- matrix(0,
    nrow = 3, ncol = 7,
    dimnames = list(
      NULL,
      c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")
    )
  )
  sim <- simulate_responses(years_old = 3, n = 3, theta = theta_mat, seed = 1)
  expect_equal(nrow(sim), 3)
})

test_that("extreme theta produces expected direction", {
  n <- 200
  factor_names <- c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")

  theta_high <- matrix(5,
    nrow = n, ncol = 7,
    dimnames = list(NULL, factor_names)
  )
  theta_low <- matrix(-5,
    nrow = n, ncol = 7,
    dimnames = list(NULL, factor_names)
  )

  sim_high <- simulate_responses(
    years_old = 3, n = n,
    theta = theta_high, seed = 1
  )
  sim_low <- simulate_responses(
    years_old = 3, n = n,
    theta = theta_low, seed = 1
  )

  item_cols <- setdiff(names(sim_high), c("id", "years_old"))
  mean_high <- mean(as.matrix(sim_high[, item_cols]))
  mean_low <- mean(as.matrix(sim_low[, item_cols]))
  expect_true(mean_high > mean_low)
})

test_that("simulate_responses errors on invalid input", {
  expect_error(simulate_responses(years_old = "three"), "numeric")
  expect_error(simulate_responses(years_old = -1), "negative")
  expect_error(simulate_responses(years_old = NA_real_), "missing")
})

test_that("simulate_responses errors on bad theta dimensions", {
  bad_theta <- matrix(0,
    nrow = 2, ncol = 7,
    dimnames = list(
      NULL,
      c("F", "GEN", "EAT", "EXT", "INT", "SLE", "SOC")
    )
  )
  expect_error(
    simulate_responses(years_old = 3, n = 5, theta = bad_theta),
    "5 rows"
  )

  bad_cols <- matrix(0, nrow = 1, ncol = 3)
  expect_error(
    simulate_responses(years_old = 3, theta = bad_cols),
    "7 columns"
  )

  no_names <- matrix(0, nrow = 1, ncol = 7)
  expect_error(
    simulate_responses(years_old = 3, theta = no_names),
    "named"
  )
})
