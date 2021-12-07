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
##' @param ci_method A method for constructing CCA confidence intervals. Needs
##'   an API-like cca_ci_asymptotic (the default)
##' @param ... additional arguments passed on to ci_method
##' @return A list of results, containing: (1) cis, a 5 dimensional array of
##'   confidence intervals, (2) truth, a 4 dimensional array that holds the true
##'   values (3) cover, a 4 dimensional array indicating which cis cover the
##'   truth, and (4) sigma, a list of length outreps that holds the generative
##'   Sigma. For each of the arrays, the first two dimensions are coordinates (p
##'   + q) and then components (K), and the last two dimensions are inreps and
##'   then outreps. For cis, the middle dimension holds lower and upper
##'   confidence bounds, respectively.
##' @export
##' @author Dan Kessler
coverage_experiment <- function(outreps, inreps, p, q, n,
                                ci_method = cca_ci_asymptotic, ...) {
  K <- min(p, q)
  coord_names <- c(paste0("X_", 1:p), paste0("Y_", 1:q))

  cis <- array(NA, c((p + q), K, 2, inreps, outreps),
    dimnames = list(
      coordinate = coord_names,
      component = 1:K,
      quantity = c("lower", "upper"),
      inreps = 1:inreps,
      outreps = 1:outreps
    )
  )

  truth <- array(NA, c((p + q), K, inreps, outreps),
    dimnames = list(
      coordinate = coord_names,
      component = 1:K,
      inreps = 1:inreps,
      outreps = 1:outreps
    )
  )

  sigma <- list()

  for (i in 1:outreps) {
    sigma[[i]] <- gen_sigma(p, q)
    fm_true <- cancor.cov(sigma[[i]], px = p)
    for (j in 1:inreps) {
      dat <- gen_data(sigma[[i]], p, q, n)
      truth[1:p, , j, i] <- fm_true$xcoef
      truth[(p + 1):(p + q), , j, i] <- fm_true$ycoef
      ci_estimates <- ci_method(dat$x, dat$y, ...)
      cis[1:p, , 1, j, i] <- ci_estimates$xcoef[, , 1]
      cis[1:p, , 2, j, i] <- ci_estimates$xcoef[, , 2]
      cis[(p + 1):(p + q), , 1, j, i] <- ci_estimates$ycoef[, , 1]
      cis[(p + 1):(p + q), , 2, j, i] <- ci_estimates$ycoef[, , 2]
    }
  }

  cover <- cis[, , 1, , ] <= drop(truth) &
    drop(truth) <= cis[, , 2, , ]

  ## Fix dims and dimnames of cover
  dim(cover) <- c((p + q), K, inreps, outreps)
  dimnames(cover) <- list(
    coordinate = coord_names,
    component = 1:K,
    inreps = 1:inreps,
    outreps = 1:outreps
  )

  return(list(cis = cis, truth = truth, cover = cover, sigma = sigma))
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
