# Tests for Weibull series distribution functions (R/wei_series.R)

test_that("dwei_series returns valid pdf values", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  # PDF should be non-negative
  expect_true(all(dwei_series(c(100, 200, 300), shapes, scales) >= 0))


  # PDF at t=0 should be 0 or finite
 expect_true(is.finite(dwei_series(0, shapes, scales)))

  # PDF for negative t should be 0
  expect_equal(dwei_series(-1, shapes, scales), 0)

  # PDF should integrate to approximately 1
  integral <- integrate(dwei_series, 0, Inf, shapes = shapes, scales = scales)$value
  expect_equal(integral, 1, tolerance = 0.01)
})

test_that("pwei_series returns valid cdf values", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  # CDF should be between 0 and 1
  expect_true(all(pwei_series(c(100, 500, 1000), shapes, scales) >= 0))
  expect_true(all(pwei_series(c(100, 500, 1000), shapes, scales) <= 1))

  # CDF at t=0 should be 0
  expect_equal(pwei_series(0, shapes, scales), 0)

  # CDF for negative t should be 0
  expect_equal(pwei_series(-1, shapes, scales), 0)

  # CDF should be monotonically increasing
  cdf_vals <- pwei_series(seq(0, 1000, 100), shapes, scales)
  expect_true(all(diff(cdf_vals) >= 0))

  # CDF approaches 1 as t -> Inf
  expect_equal(pwei_series(1e6, shapes, scales), 1, tolerance = 0.001)
})

test_that("surv_wei_series = 1 - pwei_series", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)
  t_vals <- c(100, 200, 500, 1000)

  for (t in t_vals) {
    expect_equal(
      surv_wei_series(t, shapes, scales),
      1 - pwei_series(t, shapes, scales),
      tolerance = 1e-10
    )
  }

  # Survival at t=0 should be 1
  expect_equal(surv_wei_series(0, shapes, scales), 1)

  # Survival for negative t should be 1
  expect_equal(surv_wei_series(-1, shapes, scales), 1)
})

test_that("qwei_series is inverse of pwei_series", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  # Test at various probability levels
  p_vals <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  for (p in p_vals) {
    t <- qwei_series(p, shapes, scales)
    p_recovered <- pwei_series(t, shapes, scales)
    expect_equal(p_recovered, p, tolerance = 1e-4)
  }

  # Edge cases
  expect_equal(qwei_series(0, shapes, scales), 0)
  expect_equal(qwei_series(1, shapes, scales), Inf)
})

test_that("rwei_series generates valid samples", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  set.seed(42)
  samples <- rwei_series(1000, shapes, scales)

  # All samples should be positive
  expect_true(all(samples > 0))

  # Sample mean should be close to MTTF
  mttf <- wei_series_mttf(shapes, scales)
  expect_equal(mean(samples), mttf, tolerance = mttf * 0.1)  # 10% tolerance

  # Empirical CDF should match theoretical
  emp_median <- median(samples)
  theo_median <- qwei_series(0.5, shapes, scales)
  expect_equal(emp_median, theo_median, tolerance = theo_median * 0.1)
})

test_that("hazard_wei_series is non-negative", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  # Hazard should be positive for t > 0
  expect_true(all(hazard_wei_series(c(100, 500, 1000), shapes, scales) > 0))

  # Hazard at t=0 should be NA or 0
  h0 <- hazard_wei_series(0, shapes, scales)
  expect_true(is.na(h0) || h0 == 0)
})

test_that("hazard = pdf / survival", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)
  t_vals <- c(100, 200, 500, 1000)

  for (t in t_vals) {
    h_direct <- hazard_wei_series(t, shapes, scales)
    h_formula <- dwei_series(t, shapes, scales) / surv_wei_series(t, shapes, scales)
    expect_equal(h_direct, h_formula, tolerance = 1e-10)
  }
})

test_that("wei_series_mttf returns positive value", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  mttf <- wei_series_mttf(shapes, scales)
  expect_true(mttf > 0)
  expect_true(is.finite(mttf))

  # MTTF should be less than smallest component MTTF (series system)
  # Component MTTF = scale * gamma(1 + 1/shape)
  comp_mttfs <- scales * gamma(1 + 1/shapes)
  expect_true(mttf < min(comp_mttfs))
})

test_that("wei_series_cause probabilities sum to 1", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)

  probs <- wei_series_cause(1:3, shapes = shapes, scales = scales)

  # All probabilities should be between 0 and 1
  expect_true(all(probs >= 0))
  expect_true(all(probs <= 1))

  # Sum should equal 1
  expect_equal(sum(probs), 1, tolerance = 0.01)

  # Invalid component indices should return 0
  expect_equal(wei_series_cause(0, shapes = shapes, scales = scales), 0)
  expect_equal(wei_series_cause(4, shapes = shapes, scales = scales), 0)
})

test_that("wei_series_cause with given t works", {
  shapes <- c(1.2, 1.5, 1.1)
  scales <- c(1000, 800, 900)
  t <- 500

  probs <- wei_series_cause(1:3, shapes = shapes, scales = scales, t = t)

  # All probabilities should be between 0 and 1
  expect_true(all(probs >= 0))
  expect_true(all(probs <= 1))

  # Sum should equal 1
  expect_equal(sum(probs), 1, tolerance = 1e-10)
})

test_that("distribution functions validate inputs", {
  shapes <- c(1.2, 1.5)
  scales <- c(1000, 800, 900)  # Mismatched length

  expect_error(dwei_series(100, shapes, scales))
  expect_error(pwei_series(100, shapes, scales))
  expect_error(surv_wei_series(100, shapes, scales))
  expect_error(qwei_series(0.5, shapes, scales))
  expect_error(rwei_series(10, shapes, scales))

  # Negative parameters should error
  expect_error(dwei_series(100, c(-1, 1), c(100, 100)))
  expect_error(dwei_series(100, c(1, 1), c(-100, 100)))
})
