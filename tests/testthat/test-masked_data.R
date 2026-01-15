# Tests for masked data generation (R/md_series_bernoulli_masked_component_cause.R)

test_that("md_bernoulli_cand_c1_c2_c3 assigns correct probabilities", {
  # Create test data with known component lifetimes
  df <- data.frame(
    t1 = c(100, 200, 150),  # Component 1 lifetimes
    t2 = c(200, 100, 300),  # Component 2 lifetimes
    t3 = c(300, 300, 100),  # Component 3 lifetimes
    delta = c(TRUE, TRUE, TRUE)
  )
  # Failed components: 1, 2, 3 (minimum in each row)

  p <- 0.5
  result <- md_bernoulli_cand_c1_c2_c3(df, p = p, prob = "q", comp = "t")

  # Failed component should have probability 1
  expect_equal(result$q1[1], 1)  # Row 1: comp 1 failed
  expect_equal(result$q2[2], 1)  # Row 2: comp 2 failed
  expect_equal(result$q3[3], 1)  # Row 3: comp 3 failed

  # Non-failed components should have probability p
  expect_equal(result$q2[1], p)
  expect_equal(result$q3[1], p)
  expect_equal(result$q1[2], p)
  expect_equal(result$q3[2], p)
  expect_equal(result$q1[3], p)
  expect_equal(result$q2[3], p)
})

test_that("md_bernoulli_cand_c1_c2_c3 handles right-censored data", {
  df <- data.frame(
    t1 = c(100, 200),
    t2 = c(200, 100),
    delta = c(TRUE, FALSE)  # Second observation is censored
  )

  p <- 0.5
  result <- md_bernoulli_cand_c1_c2_c3(df, p = p, prob = "q", comp = "t")

  # Censored observation should have all probabilities = 0
  expect_equal(result$q1[2], 0)
  expect_equal(result$q2[2], 0)
})

test_that("md_bernoulli_cand_c1_c2_c3 handles vector of probabilities", {
  df <- data.frame(
    t1 = c(100, 200, 300),  # Row 1: t1=100 (min), Row 2: t2=100 (min), Row 3: t3=100 (min)
    t2 = c(200, 100, 200),
    t3 = c(300, 300, 100),
    delta = c(TRUE, TRUE, TRUE)
  )

  # Different p for each observation
  p <- c(0.3, 0.5, 0.7)
  result <- md_bernoulli_cand_c1_c2_c3(df, p = p, prob = "q", comp = "t")

  # Check that failed component has probability 1
  expect_equal(result$q1[1], 1)  # Row 1: comp 1 failed
  expect_equal(result$q2[2], 1)  # Row 2: comp 2 failed
  expect_equal(result$q3[3], 1)  # Row 3: comp 3 failed

  # Check non-failed component probabilities match p
  expect_equal(result$q2[1], 0.3)  # Row 1: comp 2 not failed, p=0.3
  expect_equal(result$q3[1], 0.3)  # Row 1: comp 3 not failed, p=0.3
  expect_equal(result$q1[2], 0.5)  # Row 2: comp 1 not failed, p=0.5
  expect_equal(result$q1[3], 0.7)  # Row 3: comp 1 not failed, p=0.7
})

test_that("md_cand_sampler generates valid candidate sets", {
  # Create probability data
  df <- data.frame(
    q1 = c(1.0, 0.5, 0.5),
    q2 = c(0.5, 1.0, 0.5),
    q3 = c(0.5, 0.5, 1.0)
  )

  set.seed(42)
  result <- md_cand_sampler(df, prob = "q", candset = "x")

  # Check that candidate set columns exist
  expect_true("x1" %in% names(result))
  expect_true("x2" %in% names(result))
  expect_true("x3" %in% names(result))

  # Components with prob=1 should always be in candidate set
  expect_true(result$x1[1])
  expect_true(result$x2[2])
  expect_true(result$x3[3])

  # All values should be logical
  expect_type(result$x1, "logical")
  expect_type(result$x2, "logical")
  expect_type(result$x3, "logical")
})

test_that("md_cand_sampler respects probabilities statistically", {
  # High number of samples to test statistical properties
  n <- 1000
  df <- data.frame(
    q1 = rep(0.3, n),
    q2 = rep(0.7, n)
  )

  set.seed(42)
  result <- md_cand_sampler(df, prob = "q", candset = "x")

  # Empirical proportions should be close to true probabilities
  prop_x1 <- mean(result$x1)
  prop_x2 <- mean(result$x2)

  expect_equal(prop_x1, 0.3, tolerance = 0.05)
  expect_equal(prop_x2, 0.7, tolerance = 0.05)
})

test_that("md_bernoulli_cand_c1_c2_c3_identifiable detects issues", {
  # Good data - should return empty
  df_good <- data.frame(
    x1 = c(TRUE, TRUE, FALSE, TRUE),
    x2 = c(FALSE, TRUE, TRUE, TRUE),
    x3 = c(TRUE, FALSE, TRUE, FALSE)
  )
  reasons <- md_bernoulli_cand_c1_c2_c3_identifiable(df_good, candset = "x")
  expect_equal(length(reasons), 0)

  # Bad: component never in candidate set
  df_never <- data.frame(
    x1 = c(FALSE, FALSE, FALSE),
    x2 = c(TRUE, TRUE, TRUE)
  )
  reasons <- md_bernoulli_cand_c1_c2_c3_identifiable(df_never, candset = "x")
  expect_true(length(reasons) > 0)

  # Bad: component always in candidate set
  df_always <- data.frame(
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, TRUE)
  )
  reasons <- md_bernoulli_cand_c1_c2_c3_identifiable(df_always, candset = "x")
  expect_true(any(grepl("always", reasons, ignore.case = TRUE)))

  # Bad: all candidate sets have all components
  df_all <- data.frame(
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, TRUE, TRUE)
  )
  reasons <- md_bernoulli_cand_c1_c2_c3_identifiable(df_all, candset = "x")
  expect_true(length(reasons) > 0)
})

test_that("md_decorate_component_cause_k adds correct cause", {
  df <- data.frame(
    t1 = c(100, 200, 150),
    t2 = c(200, 100, 300),
    t3 = c(300, 300, 100)
  )

  result <- md_decorate_component_cause_k(df, comp = "t")

  expect_true("k" %in% names(result))
  expect_equal(result$k, c(1, 2, 3))  # Minimum components
})

test_that("md_cand_sampler validates input", {
  # Non-data frame should error
  expect_error(md_cand_sampler(list(q1 = c(0.5, 0.5)), prob = "q", candset = "x"))

  # Missing probability columns should error
  df <- data.frame(a = c(1, 2), b = c(3, 4))
  expect_error(md_cand_sampler(df, prob = "q", candset = "x"))

  # Empty data frame should return empty
  df_empty <- data.frame(q1 = numeric(0), q2 = numeric(0))
  result <- md_cand_sampler(df_empty, prob = "q", candset = "x")
  expect_equal(nrow(result), 0)
})

test_that("md_bernoulli_cand_c1_c2_c3 validates input", {
  # Missing component columns should error
  df <- data.frame(
    a = c(100, 200),
    b = c(200, 100),
    delta = c(TRUE, TRUE)
  )
  expect_error(md_bernoulli_cand_c1_c2_c3(df, p = 0.5, comp = "t"))

  # Missing right censoring column should error
  df <- data.frame(
    t1 = c(100, 200),
    t2 = c(200, 100)
  )
  expect_error(md_bernoulli_cand_c1_c2_c3(
    df, p = 0.5, comp = "t",
    right_censoring_indicator = "delta"
  ))

  # Empty data frame should return empty
  df_empty <- data.frame(t1 = numeric(0), t2 = numeric(0), delta = logical(0))
  result <- md_bernoulli_cand_c1_c2_c3(df_empty, p = 0.5, comp = "t")
  expect_equal(nrow(result), 0)
})

test_that("full masking pipeline works end-to-end", {
  # Generate component lifetimes
  set.seed(42)
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)
  n <- 50

  # Create component lifetimes
  t1 <- stats::rweibull(n, shape = shapes[1], scale = scales[1])
  t2 <- stats::rweibull(n, shape = shapes[2], scale = scales[2])
  t3 <- stats::rweibull(n, shape = shapes[3], scale = scales[3])

  df <- data.frame(t1 = t1, t2 = t2, t3 = t3, delta = rep(TRUE, n))

  # Add masking probabilities
  df <- md_bernoulli_cand_c1_c2_c3(df, p = 0.5, prob = "q", comp = "t")

  # Generate candidate sets
  df <- md_cand_sampler(df, prob = "q", candset = "x")

  # Add system lifetime (minimum of components)
  df$t <- pmin(df$t1, df$t2, df$t3)

  # Run MLE
  theta0 <- c(1, 1000, 1, 1000, 1, 1000)
  result <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 500L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )

  # MLE should converge
  expect_equal(result$convergence, 0)

  # Parameters should be in reasonable range
  expect_true(all(result$par > 0))
  expect_true(all(result$par[seq(1, 6, 2)] < 5))  # Shapes reasonable
  expect_true(all(result$par[seq(2, 6, 2)] < 5000))  # Scales reasonable
})
