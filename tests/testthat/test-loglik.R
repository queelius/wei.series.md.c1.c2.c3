# Tests for likelihood functions (R/md_wei_series_c1_c2_c3.R)

test_that("loglik_wei_series_md_c1_c2_c3 returns finite value for valid data", {
  # Use bundled guo data
  df <- guo_weibull_series_md$data
  theta <- guo_weibull_series_md$mle

  ll <- loglik_wei_series_md_c1_c2_c3(df, theta)

  expect_true(is.finite(ll))
  expect_true(ll < 0)  # Log-likelihood should be negative
})

test_that("loglik matches Guo et al. published value", {
  df <- guo_weibull_series_md$data
  theta <- guo_weibull_series_md$mle

  ll <- loglik_wei_series_md_c1_c2_c3(df, theta)

  # Should match the published log-likelihood
  expect_equal(ll, guo_weibull_series_md$loglike, tolerance = 0.01)
})

test_that("loglik is maximized at MLE", {
  df <- guo_weibull_series_md$data
  theta_mle <- guo_weibull_series_md$mle
  ll_mle <- loglik_wei_series_md_c1_c2_c3(df, theta_mle)

  # Perturb parameters and check that log-likelihood decreases
  perturbations <- list(
    c(1.5, 1000, 1.2, 900, 1.1, 850),  # Different shapes
    c(1.2, 800, 1.2, 800, 1.2, 800),   # Different scales
    c(1.0, 1000, 1.0, 1000, 1.0, 1000) # Uniform parameters
  )

  for (theta_perturbed in perturbations) {
    ll_perturbed <- loglik_wei_series_md_c1_c2_c3(df, theta_perturbed)
    expect_true(ll_mle >= ll_perturbed)
  }
})

test_that("score function is approximately zero at MLE", {
  df <- guo_weibull_series_md$data
  theta_mle <- guo_weibull_series_md$mle

  score <- score_wei_series_md_c1_c2_c3(df, theta_mle)

  # Score should be close to zero at MLE
  expect_true(all(abs(score) < 1))  # Relaxed tolerance for numerical MLE
})

test_that("score_slow matches score (vectorized version)", {
  df <- guo_weibull_series_md$data
  theta <- c(1.2, 900, 1.1, 850, 1.15, 800)

  score_fast <- score_wei_series_md_c1_c2_c3(df, theta)
  score_slow <- score_wei_series_md_c1_c2_c3_slow(df, theta)

  expect_equal(score_fast, score_slow, tolerance = 1e-10)
})

test_that("hessian is negative definite at MLE", {
  df <- guo_weibull_series_md$data
  theta_mle <- guo_weibull_series_md$mle

  H <- hessian_wei_series_md_c1_c2_c3(df, theta_mle)

  # Hessian should be symmetric
  expect_equal(H, t(H), tolerance = 1e-6)

  # All eigenvalues should be negative (negative definite)
  eigenvalues <- eigen(H, symmetric = TRUE)$values
  expect_true(all(eigenvalues < 0))
})

test_that("loglik handles right-censored data correctly", {
  # Create simple test data with some censored observations
  df <- data.frame(
    t = c(100, 200, 300, 400, 500),
    x1 = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    x2 = c(FALSE, TRUE, TRUE, FALSE, TRUE),
    delta = c(TRUE, TRUE, TRUE, FALSE, FALSE)  # Last two are censored
  )

  theta <- c(1.5, 500, 1.2, 400)

  ll <- loglik_wei_series_md_c1_c2_c3(df, theta)
  expect_true(is.finite(ll))

  # Without censoring indicator, should assume all observed
  df_no_delta <- df[, c("t", "x1", "x2")]
  ll_no_censor <- loglik_wei_series_md_c1_c2_c3(
    df_no_delta, theta,
    right_censoring_indicator = "nonexistent"
  )
  expect_true(is.finite(ll_no_censor))

  # Log-likelihood with censoring should be different
  expect_false(ll == ll_no_censor)
})

test_that("loglik validates input parameters", {
  df <- guo_weibull_series_md$data

  # Wrong theta length should error
  expect_error(loglik_wei_series_md_c1_c2_c3(df, c(1, 1, 1, 1)))  # Too short

  # Missing lifetime column should error
  df_bad <- df
  names(df_bad)[names(df_bad) == "t"] <- "time"
  expect_error(loglik_wei_series_md_c1_c2_c3(df_bad, guo_weibull_series_md$mle))

  # Empty data frame should error
  expect_error(loglik_wei_series_md_c1_c2_c3(df[0, ], guo_weibull_series_md$mle))
})

test_that("score validates input parameters", {
  df <- guo_weibull_series_md$data

  # Missing lifetime column should error
  df_bad <- df
  names(df_bad)[names(df_bad) == "t"] <- "time"
  expect_error(score_wei_series_md_c1_c2_c3(df_bad, guo_weibull_series_md$mle))

  # Empty data frame should error
  expect_error(score_wei_series_md_c1_c2_c3(df[0, ], guo_weibull_series_md$mle))
})

test_that("loglik with custom column names works", {
  df <- guo_weibull_series_md$data

  # Rename columns
  df_renamed <- df
  names(df_renamed)[names(df_renamed) == "t"] <- "lifetime"
  names(df_renamed)[names(df_renamed) == "delta"] <- "observed"
  names(df_renamed) <- gsub("^x", "cand", names(df_renamed))

  theta <- guo_weibull_series_md$mle

  ll <- loglik_wei_series_md_c1_c2_c3(
    df_renamed, theta,
    lifetime = "lifetime",
    candset = "cand",
    right_censoring_indicator = "observed"
  )

  ll_original <- loglik_wei_series_md_c1_c2_c3(df, theta)

  expect_equal(ll, ll_original)
})
