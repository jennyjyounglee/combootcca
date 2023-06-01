
test_that("cancor_cov is a proper inverse of cancor_inv_cov regardless of dimensions", {
  p <- 2
  q <- 5
  K <- min(p, q)

  fm1 <- list()
  fm1$cor <- sort(runif(K), decreasing = TRUE)

  fm1$xcoef <- matrix(rnorm(p * K), p, K)
  fm1$ycoef <- matrix(rnorm(q * K), q, K)

  fm1 <- cca_align_posdiag(fm1)

  Sigma <- do.call(cancor_inv_cov, fm1)

  fm2 <- cancor_cov(Sigma, p, align = cca_align_posdiag)

  expect_equal(fm1, fm2)
})


test_that("cancor_inv_cov is an inverse of cancor_inv_cov WHEN p = q", {
  p <- 5
  q <- 5
  K <- min(p, q)

  Q <- randortho_fixed(p + q)
  eigs <- sort(runif(p + q), decreasing = FALSE)

  Sigma <- Q %*% diag(eigs) %*% t(Q)

  fm1 <- cancor_cov(Sigma, p)

  Sigma2 <- do.call(cancor_inv_cov, fm1)

  expect_equal(Sigma, Sigma2)
})

