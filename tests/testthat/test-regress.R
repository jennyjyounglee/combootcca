test_that("lm_scalefix gives unit variance predictions", {
  n <- 100
  p <- 5

  y <- rnorm(n)
  x <- matrix(rnorm(n * p), n, p)

  fm <- lm_scalefix(lm(y ~ x))

  expect_equal(var(predict(fm)), 1)

})

test_that("CCA and Regression Agree (when appropriately rescaled)", {

  n <- 100
  p <- 5

  y <- rnorm(n)
  x <- matrix(rnorm(n * p), n, p)

  fm_regress <- lm_scalefix(lm(y ~ x))
  fm_cca <- cancor_scaled(x, y)

  sign_flip <- sign(fm_cca$xcoef[1, 1]) != sign(coef(fm_regress)[2])
  if (sign_flip) {
    fm_cca$xcoef[, 1] <- -1 * fm_cca$xcoef[, 1]
  }

  expect_equal(fm_cca$xcoef[, 1], unname(coef(fm_regress)[-1]))

})
