# Tests for right censoring functions (R/md_series_lifetime_right_censoring.R)

test_that("md_series_lifetime_right_censoring computes system lifetime", {
  df <- data.frame(
    t1 = c(100, 200, 150),
    t2 = c(200, 100, 300),
    t3 = c(300, 300, 100)
  )

  result <- md_series_lifetime_right_censoring(df, tau = Inf, comp = "t")

  # System lifetime is minimum of components
  expect_equal(result$t, c(100, 100, 100))

  # All should be observed (no censoring when tau = Inf)
  expect_true(all(result$delta))
})

test_that("md_series_lifetime_right_censoring handles censoring correctly", {
  df <- data.frame(
    t1 = c(100, 200, 400),
    t2 = c(200, 100, 500)
  )

  # Right censor at tau = 150
  result <- md_series_lifetime_right_censoring(df, tau = 150, comp = "t")

  # Row 1: min(100,200) = 100 < 150, observed
  expect_equal(result$t[1], 100)
  expect_true(result$delta[1])

  # Row 2: min(200,100) = 100 < 150, observed
  expect_equal(result$t[2], 100)
  expect_true(result$delta[2])

  # Row 3: min(400,500) = 400 > 150, censored at tau
  expect_equal(result$t[3], 150)
  expect_false(result$delta[3])
})

test_that("md_series_lifetime_right_censoring validates input", {
  df <- data.frame(a = c(1, 2), b = c(3, 4))

  # Missing component columns should error
  expect_error(md_series_lifetime_right_censoring(df, tau = 100, comp = "t"))
})

test_that("md_series_lifetime_right_censoring handles vector tau", {
  df <- data.frame(
    t1 = c(100, 200, 300),
    t2 = c(200, 100, 400)
  )

  # Different tau for each observation
  result <- md_series_lifetime_right_censoring(df, tau = c(50, 150, 500), comp = "t")

  # Row 1: min(100,200) = 100 > 50, censored at 50
  expect_equal(result$t[1], 50)
  expect_false(result$delta[1])

  # Row 2: min(200,100) = 100 < 150, observed
  expect_equal(result$t[2], 100)
  expect_true(result$delta[2])

  # Row 3: min(300,400) = 300 < 500, observed
  expect_equal(result$t[3], 300)
  expect_true(result$delta[3])
})
