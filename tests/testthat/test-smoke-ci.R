test_that("All methods for CCA confidence intervals run without errors", {
  n <- 100
  p <- 3
  q <- 4

  x <- matrix(rnorm(n * p), n, p)
  y <- matrix(rnorm(n * q), n, q)

  expect_error(cca_ci_asymptotic(x, y), NA)
  expect_error(cca_ci_regression(x, y), NA)
  expect_error(cca_ci_bootstrap_abs(x, y, parametric = FALSE), NA)
  expect_error(cca_ci_bootstrap_abs(x, y, parametric = TRUE), NA)
})
})
