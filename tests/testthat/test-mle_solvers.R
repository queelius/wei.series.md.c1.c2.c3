# Tests for MLE solvers (R/md_mle_wei_series_c1_c2_c3.R)

test_that("mle_lbfgsb recovers known MLE", {
  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  result <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 1000L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )

  # Check convergence
  expect_equal(result$convergence, 0)

  # Check that MLE matches Guo et al.
  expect_equal(result$par, guo_weibull_series_md$mle, tolerance = 0.01)

  # Check that log-likelihood matches
  expect_equal(-result$value, -guo_weibull_series_md$loglike, tolerance = 0.01)
})

test_that("mle_lbfgsb returns hessian when requested",
{
  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  result_with_hess <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    hessian = TRUE,
    control = list(maxit = 500L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )

  result_without_hess <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    hessian = FALSE,
    control = list(maxit = 500L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )

  expect_true(!is.null(result_with_hess$hessian))
  expect_true(is.null(result_without_hess$hessian))
})

test_that("mle_nelder finds approximate MLE", {
  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  # Suppress NaN warnings (expected when optimizer explores invalid parameter space)
  result <- suppressWarnings(mle_nelder_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 2000L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  ))

  # Nelder-Mead may not be as accurate, use relaxed tolerance
  expect_equal(result$par, guo_weibull_series_md$mle, tolerance = 0.1)
})

test_that("mle_newton finds approximate MLE", {
  df <- guo_weibull_series_md$data
  # Start closer to the true MLE for Newton's method stability
  theta0 <- c(1.2, 900, 1.1, 850, 1.1, 800)

  result <- mle_newton_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 200L, lr = 0.5)
  )

  # Newton may not converge perfectly, check it improves log-likelihood
  ll_initial <- loglik_wei_series_md_c1_c2_c3(df, theta0)
  expect_true(result$value >= ll_initial)
})

test_that("mle_grad improves log-likelihood", {
  df <- guo_weibull_series_md$data
  theta0 <- c(1.2, 900, 1.1, 850, 1.1, 800)

  result <- mle_grad_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 100L, lr = 0.01)
  )

  ll_initial <- loglik_wei_series_md_c1_c2_c3(df, theta0)
  expect_true(result$value >= ll_initial)
})

test_that("mle_sann runs without error", {
  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  # SANN is stochastic and slow, just test it runs
  # Suppress NaN warnings (expected when optimizer explores invalid parameter space)
  result <- suppressWarnings(mle_sann_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 100L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  ))

  expect_true(!is.null(result$par))
  expect_true(is.finite(result$value))
})

test_that("all solvers respect parameter bounds", {
  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  # L-BFGS-B explicitly enforces bounds
  result_lbfgsb <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 500L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )
  expect_true(all(result_lbfgsb$par > 0))

  # Nelder-Mead (suppress NaN warnings)
  result_nelder <- suppressWarnings(mle_nelder_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 500L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  ))
  expect_true(all(result_nelder$par > 0))
})

test_that("mle_lbfgsb works with algebraic.mle wrapper", {
  skip_if_not_installed("algebraic.mle")

  df <- guo_weibull_series_md$data
  theta0 <- c(1, 1, 1, 1, 1, 1)

  result <- mle_lbfgsb_wei_series_md_c1_c2_c3(
    df = df,
    theta0 = theta0,
    control = list(maxit = 1000L, parscale = c(1, 1000, 1, 1000, 1, 1000))
  )

  # Wrap with mle_numerical
  fit <- algebraic.mle::mle_numerical(result)

  # Check that algebraic.mle methods work
  expect_equal(length(algebraic.mle::params(fit)), 6)
  expect_true(is.finite(algebraic.mle::loglik_val(fit)))

  # Confidence intervals should work
  ci <- confint(fit)
  expect_equal(nrow(ci), 6)
  expect_equal(ncol(ci), 2)

  # All lower bounds should be less than upper bounds
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("solvers handle different starting points", {
  df <- guo_weibull_series_md$data

  # Multiple starting points
  starts <- list(
    c(0.5, 500, 0.5, 500, 0.5, 500),
    c(2, 1500, 2, 1500, 2, 1500),
    c(1, 1000, 1, 1000, 1, 1000)
  )

  results <- lapply(starts, function(theta0) {
    mle_lbfgsb_wei_series_md_c1_c2_c3(
      df = df,
      theta0 = theta0,
      control = list(maxit = 1000L, parscale = c(1, 1000, 1, 1000, 1, 1000))
    )
  })

  # All should converge to similar MLE
  for (i in 2:length(results)) {
    expect_equal(results[[1]]$par, results[[i]]$par, tolerance = 0.1)
  }
})
