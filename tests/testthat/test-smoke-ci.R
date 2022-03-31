test_that("All methods for CCA confidence intervals run without errors", {
  n <- 100 # needs to be less than nboots for BCA to work
  p <- 3
  q <- 4

  x <- matrix(rnorm(n * p), n, p)
  y <- matrix(rnorm(n * q), n, q)

  expect_error(cca_ci_asymptotic(x, y), NA)
  expect_error(cca_ci_regression(x, y), NA)
  expect_error(cca_ci_absboot(x, y, parametric = FALSE), NA)
  expect_error(cca_ci_absboot(x, y, parametric = TRUE), NA)
  suppressWarnings(expect_error(cca_ci_boot(x, y, boot_type = c("norm", "basic", "perc")), NA))
  expect_error(cca_ci_boot(x, y, boot_type = c("norm", "basic", "perc")), NA) # use a subset of boot_types
  suppressWarnings(expect_error(cca_ci_boot(x, y, ncpus = 5, boot_type = c("norm", "basic", "perc")), NA)) # test multicore support
})

test_that("Standard Metrics Work", {

  n <- 1e6
  p <- 3
  q <- 3

  Sigma <- gen_sigma(p, q)
  fm_true <- cancor_cov(Sigma, px = p)

  dat <- gen_data(Sigma, p, q, n)

  fm_hat <- cancor_scaled(dat$x, dat$y)

  cis <- cca_ci_asymptotic(dat$x, dat$y, level = 0.95)

  expect_error(cca_metric_standard(fm_true, cis), NA)


})

