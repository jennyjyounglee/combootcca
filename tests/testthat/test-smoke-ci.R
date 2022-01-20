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
  suppressWarnings(expect_error(cca_ci_boot(x, y), NA))
})

test_that("Coverage Assessments Works", {

  n <- 1e3
  p <- 20
  q <- 20

  Sigma <- gen_sigma(p, q)
  fm_true <- cancor_cov(Sigma, px = p)

  dat <- gen_data(Sigma, p, q, n)

  fm_hat <- cancor_scaled(dat$x, dat$y)

  cis <- cca_ci_regression(dat$x, dat$y, align = cca_align_greedy_cosx, ref = fm_hat)

  cca_ci_coverage_possibilities(fm_true, cis)

  cca_ci_coverage_pangloss(fm_true, cis)

  cca_ci_coverage_signflip(fm_true, cis)


})
