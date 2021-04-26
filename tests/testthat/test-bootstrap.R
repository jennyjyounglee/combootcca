test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("CCA on Covariance and Data Agree", {
    n <- 100
    px <- 5
    py <- 6
    Sigma <- generate.sigma(px, py)
    dat <- generate.data(n, Sigma, px)
    X <- dat$X
    Y <- dat$Y
    fm1 <- cancor(dat$X, dat$Y)
    Sigma.hat <- cov(do.call(cbind, dat))
    fm2 <- cancor.cov(Sigma.hat, px)

    ## check constraint
    cov(X %*% fm1$xcoef)
    cov(Y %*% fm1$ycoef)

    Xc <- scale(X, center=TRUE, scale=FALSE)
    Yc <- scale(Y, center=TRUE, scale=FALSE)

    t(Xc %*% fm1$xcoef) %*% (Xc %*% fm1$xcoef)
    
    cov(X %*% fm2$xcoef)
    cov(Y %*% fm2$ycoef)

})
