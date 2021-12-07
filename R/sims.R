## code to support simulation studies of CCA

sim.coverage.check <- function(n, px, py) {
  Sigma <- generate.sigma(px, py)
  dat <- generate.data(n, Sigma, px)
  X <- dat$X
  Y <- dat$Y

  pop.fm <- cancor.cov(Sigma, px)

  boot.cis <- bootstrapcca(X, Y) }
##' @title Conduct asymptotic coverage experiment
##' @param outreps Each "outer replication" draws a new value of sigma
##' @param inreps Each "inner replication" is a repetition with new data for the
##'   same value of sigma
##' @param p The dimension of X
##' @param q The dimension of Y
##' @param n How many datapoints to draw in each sample
##' @return A big array with all of the results. The array has the following
##'   dimensions: outreps by inreps by (p + q) by K by 3. The last dimension
##'   holds the lower CI estimates, upper CI estimates, and true values,
##'   respectively.
##' @export
##' @author Dan Kessler
coverage_experiment <- function(outreps, inreps, p, q, n) {
  K <- min(p, q)
  output <- array(NA, c(outreps, inreps, (p + q), K, 3),
    dimnames = list(
      outreps = NULL,
      inreps = NULL,
      coordinate = NULL,
      component = NULL,
      quantity = c("lower", "upper", "truth")
    )
  )

  for (i in 1:outreps) {
    sigma <- gen_sigma(p, q)
    fm_true <- cancor.cov(sigma, px = p)
    for (j in 1:inreps) {
      dat <- gen_data(sigma, p, q, n)
      output[i, j, 1:p, , 3] <- fm_true$xcoef
      output[i, j, (p + 1):(p + q), , 3] <- fm_true$ycoef
      ci_asymptotic <- cca_ci_asymptotic(dat$x, dat$y)
      output[i, j, 1:p, , 1] <- ci_asymptotic$xcoef[, , 1]
      output[i, j, 1:p, , 2] <- ci_asymptotic$xcoef[, , 2]
      output[i, j, (p + 1):(p + q), , 1] <- ci_asymptotic$ycoef[, , 1]
      output[i, j, (p + 1):(p + q), , 2] <- ci_asymptotic$ycoef[, , 2]
    }
  }
  return(output)
}

coverage_asymptotic <- function(n, p, q, reps) {
  sigma <- gen_sigma(p, q)

  sigma_r <- chol(sigma)
  fm_true <- cancor.cov(sigma, px = p)

  cover <- list()
  cover$x <- array(NA, c(dim(fm_true$xcoef), reps))
  cover$y <- array(NA, c(dim(fm_true$ycoef), reps))

  for (i in 1:reps) {
    sigma <- gen_sigma(p, q)
    sigma_r <- chol(sigma)
    fm_true <- cancor.cov(sigma, px = p)

    newdata <- t(t(sigma_r) %*% matrix(rnorm(n * (p + q)), p + q, n))
    x <- newdata[, 1:p]
    y <- newdata[, (p + 1):(p + q)]

    ci_asymptotic <- cca_ci_asymptotic(x, y)

    cover$x[, , i] <- ci_asymptotic$xcoef[, , 1] < fm_true$xcoef &
      fm_true$xcoef < ci_asymptotic$xcoef[, , 2]

    cover$y[, , i] <- ci_asymptotic$ycoef[, , 1] < fm_true$ycoef &
      fm_true$ycoef < ci_asymptotic$ycoef[, , 2]
  }
  return(cover)
}

gen_data <- function(Sigma, p, q, n) {
  Sigma_r <- chol(Sigma)

  newdata <- t(t(Sigma_r) %*% matrix(rnorm(n * (p + q)), p + q, n))
  x <- newdata[, 1:p]
  y <- newdata[, (p + 1):(p + q)]
  res <- list(x = x, y = y)
  return(res)
}

gen_sigma <- function(p, q) {
  sxx <- diag(p)
  syy <- diag(q)

  sxx_sqrt <- expm::sqrtm(sxx)
  syy_sqrt <- expm::sqrtm(syy)

  qx <- randortho_fixed(p)
  qy <- randortho_fixed(q)

  rho <- sort(runif(min(p, q)), decreasing = TRUE)
  s <- matrix(0, q, p)
  diag(s) <- rho

  syx <- syy_sqrt %*% qy %*% s %*% t(qx) %*% sxx_sqrt

  sigma <- matrix(0, p + q, p + q)
  sigma[1:p, 1:p] <- sxx
  sigma[(p + 1):(p + q), (p + 1):(p + q)] <- syy
  sigma[(p + 1):(p + q), 1:p] <- syx
  sigma[1:p, (p + 1):(p + q)] <- t(syx)

  return(sigma)
}

## like pracma::randortho but without the bug
randortho_fixed <- function(n, type = c("orthonormal", "unitary")) {
    stopifnot(is.numeric(n), length(n) == 1,
              floor(n) == ceiling(n), n >= 1)
    if (n == 1)
        return(matrix(1, 1, 1))

    type <- match.arg(type)
    if (type == "orthonormal") {
        z <- pracma::randn(n, n) / sqrt(2.0)
    } else {
        z <- (pracma::randn(n, n) + 1i * pracma::randn(n, n)) / sqrt(2.0)
    }

    # QR decomposition for real or complex matrices
    Z <- qr(z)
    q <- qr.Q(Z)
    r <- qr.R(Z)

    d <- diag(r)
    ph <- d / abs(d)
    q %*% diag(ph)
}
