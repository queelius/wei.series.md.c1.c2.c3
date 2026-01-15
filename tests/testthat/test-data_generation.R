# Tests for data generation functions

test_that("generate_guo_weibull_table_2_data creates valid data", {
  set.seed(42)
  df <- generate_guo_weibull_table_2_data(n = 10)

  # Should have correct structure
  expect_true("t" %in% names(df))
  expect_true("delta" %in% names(df))
  expect_true("x1" %in% names(df))
  expect_true("x2" %in% names(df))
  expect_true("x3" %in% names(df))

  # Should have correct number of rows
  expect_equal(nrow(df), 10)

  # Lifetimes should be positive
  expect_true(all(df$t > 0))

  # Candidate sets should be logical
  expect_type(df$x1, "logical")
  expect_type(df$x2, "logical")
  expect_type(df$x3, "logical")

  # At least one component should be in each candidate set for observed failures
  observed <- df[df$delta, ]
  cand_sums <- rowSums(observed[, c("x1", "x2", "x3")])
  expect_true(all(cand_sums >= 1))
})

test_that("generate_guo_weibull_table_2_data handles custom parameters", {
  set.seed(42)
  df <- generate_guo_weibull_table_2_data(
    n = 5,
    shapes = c(1, 1, 1),
    scales = c(500, 500, 500),
    p = 0.8,
    tau = 1000
  )

  expect_equal(nrow(df), 5)
  expect_true(all(df$t <= 1000))  # Respect censoring time
})

test_that("generated data can be used for MLE", {
  set.seed(42)
  df <- generate_guo_weibull_table_2_data(n = 50)

  # Should be able to compute log-likelihood
  theta <- c(1.2, 900, 1.1, 850, 1.1, 800)
  ll <- loglik_wei_series_md_c1_c2_c3(df, theta)

  expect_true(is.finite(ll))
  expect_true(ll < 0)
})
