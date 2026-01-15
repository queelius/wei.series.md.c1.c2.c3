# Tests for Weibull component functions (R/wei_series_comp.R)

test_that("hazard_wei returns correct values", {
  # Weibull hazard: h(t) = (shape/scale) * (t/scale)^(shape-1)
  shape <- 2
  scale <- 100
  t <- 50

  expected <- shape / scale * (t / scale)^(shape - 1)
  actual <- hazard_wei(t, shape, scale)

  expect_equal(actual, expected)
})

test_that("hazard_wei is vectorized over t", {
  shape <- 1.5
  scale <- 100
  t_vals <- c(10, 50, 100)

  results <- hazard_wei(t_vals, shape, scale)

  expect_equal(length(results), 3)
  expect_true(all(results > 0))
})

test_that("wei_mttf returns correct value", {
  # Weibull MTTF = scale * gamma(1 + 1/shape)
  shape <- 2
  scale <- 100

  expected <- scale * gamma(1 + 1/shape)
  actual <- wei_mttf(shape, scale)

  expect_equal(actual, expected)
})

test_that("wei_mttf equals scale when shape = 1 (exponential)", {
  # For exponential (shape=1), MTTF = scale * gamma(2) = scale * 1 = scale
  shape <- 1
  scale <- 100

  expect_equal(wei_mttf(shape, scale), scale)
})
