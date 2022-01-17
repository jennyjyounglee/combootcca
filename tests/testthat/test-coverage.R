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
