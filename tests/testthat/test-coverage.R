test_that("Coverage is computed correctly", {
  truth <- c(1, 2, 3, 4)
  cis <- matrix(c(
    0, 2,
    1, 3,
    0, 1,
    3.5, 4.1
  ),
  4, 2,
  byrow = TRUE
  )

  expect_equal(cca_ci_coverage1(truth, cis), 3 / 4)
})

test_that("Signal detection low level functions work", {
  true <- diag(2)
  ci_lower <- matrix(c(-1, -1, 1, 1), 2, 2)
  ci_upper <- ci_lower + 2
  ci <- abind::abind(ci_lower, ci_upper, along = 3)

  out <- cca_ci_sigdet1(true, ci)

  expect_equal(out$value, c(1, 0, 2, 3))
})
