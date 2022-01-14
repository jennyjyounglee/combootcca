test_that("All methods for CCA confidence intervals run without errors", {
  n <- 100
  p <- 3
  q <- 4

  x <- matrix(rnorm(n * p), n, p)
  y <- matrix(rnorm(n * q), n, q)

  expect_error(cca_ci_asymptotic(x, y), NA)
  expect_error(cca_ci_regression(x, y), NA)
  expect_error(cca_ci_absboot(x, y, parametric = FALSE), NA)
  expect_error(cca_ci_absboot(x, y, parametric = TRUE), NA)
})

test_that("Inner BOOT Function Works", {

 problem <- bt_problem_std_fun(p = 2, q = 3, n = 1000, inreps = 3L)

 expect_error(bt_algo_boot_inner(problem$inreps[[1]], 1000,
   align = cca_align_posdiag,
   level = .90
 ), NA)

})
