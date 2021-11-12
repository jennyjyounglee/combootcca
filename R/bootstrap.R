##' Bootstrap data and
##'
##' .. content for \details{} ..
##' @title bootstrap CCA with Sign Fixing
##' @param X Data matrix of size N by p
##' @param Y Data matrix of size N by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param nboots Number of bootstrap sample to draw
##' @return
##' @author Daniel Kessler
##' @export
cca_ci_bootstrap <- function(x, y, level = .05, nboots = 1e3) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  K <- min(p, q)

  fm <- cancor_scalefix(cancor(x, y), n)
  K <- length(fm$cor)
  
  rho_hat <- fm$cor
  xcoef_hat <- fm$xcoef[, 1:K]
  ycoef_hat <- fm$ycoef[, 1:K]

  rho_boot <- array(dim = c(length(rho_hat), nboots))
  xcoef_boot <- array(dim = c(dim(xcoef_hat), nboots))
  ycoef_boot <- array(dim = c(dim(ycoef_hat), nboots))

  for (i in 1:nboots) {
    idx <- sample(n, replace = TRUE)
    fm_boot <- cancor_scalefix(cancor(x[idx, ], y[idx, ]), n)
    rho_boot[, i] <- fm_boot$cor
    xcoef_boot[, , i] <- fm_boot$xcoef[, 1:K]
    ycoef_boot[, , i] <- fm_boot$ycoef[, 1:K]
  }

  rho_dist <- abs(sweep(rho_boot, c(1), rho_hat))
  xcoef_dist <- abs(sweep(xcoef_boot, c(1, 2), xcoef_hat))
  ycoef_dist <- abs(sweep(ycoef_boot, c(1, 2), ycoef_hat))

  rho_t <- apply(rho_dist, 1, quantile, probs = level)
  xcoef_t <- apply(xcoef_dist, c(1, 2), quantile, probs = level)
  ycoef_t <- apply(ycoef_dist, c(1, 2), quantile, probs = level)

  rho_ci <- array(c(rho_hat - rho_t, rho_hat, rho_hat + rho_t),
    dim = c(K, 3),
    dimnames = list(Component = 1:K, Interval = c("Lower", "Estimate", "Upper"))
  )

  xcoef_ci <- array(NA, c(p, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
  )

  ycoef_ci <- array(NA, c(q, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
  )

  xcoef_ci[, , 1] <- xcoef_hat - xcoef_t
  xcoef_ci[, , 1] <- xcoef_hat + xcoef_t

  ycoef_ci[, , 1] <- ycoef_hat - ycoef_t
  ycoef_ci[, , 1] <- ycoef_hat + ycoef_t


  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}

##' Correct sign ambiguity in canonical correlation analysis
##'
##'
##' Canonical Correlation Analysis suffers from sign ambiguity in the
##' estimated coefficients. This is because Corr(U, V) = Corr(-U, -V).
##' This function adopts the convention that each canonical
##' coefficient pair (xcoef, ycoef) should satisfy the condition that
##' its maximal element (in absolute value) has positive sign. You may
##' wish to standardize the columns of X and Y prior to initially
##' fitting CCA so that differences in variable scaling do not
##' dominate the identity of the maximal element.
##'
##' Note: This only cares about and tries to fix sign ambiguities in
##' the first K=min(px,py) components.
##' @title Fix Cancor Signs
##' @param fm A fitted object returned by cancor
##' @return Same object as returned by cancor after sign-flipping per
##'     the identifiability condition discussed in Details.
##' @author Daniel Kessler
cancor_signfix <- function(fm) {
  K <- min(ncol(fm$xcoef), ncol(fm$ycoef))

  theta <- rbind(fm$xcoef[, 1:K], fm$ycoef[, 1:K])

  maxes <- apply(theta, 2, absmax)

  flip <- diag(sign(maxes))

  fm$xcoef[, 1:K] <- fm$xcoef[, 1:K] %*% flip
  fm$ycoef[, 1:K] <- fm$ycoef[, 1:K] %*% flip

  return(fm)
}

##' Adjust canonical correlation analysis results so that canonical
##' variates have unit variance
##'
##' Cancor uses the Golub and Van Loan algorithm, which yields
##' canonical variates that have unit norm rather than unit variance.
##' As a result, the canonical variates will have (empirical) variance
##' of 1/(N-1), where N is the number of observations. This is
##' undesirable if we wish to perform inference on the weights in
##' xcoef or ycoef, because their scale will vary with N.
##'
##' This function simply multiplies both xcoef and ycoef by sqrt(N-1)
##' so that the resulting canonical variates will indeed have unit variance.
##' @title Fix Cancor Scaling
##' @param fm A fitted object returned by cancor
##' @param N The number of observations
##' @return A modified version of fm
##' @author Daniel Kessler
cancor_scalefix <- function(fm, N) {
  fm$xcoef <- sqrt(N - 1) * fm$xcoef
  fm$ycoef <- sqrt(N - 1) * fm$ycoef
  return(fm)
}

##' Compute the Canonical Correlation Analysis solution based on a
##' covariance matrix
##'
##' The Canonical Correlation Analysis problem can be solved in a
##' variety of ways. The approach taken by the cancor function in R
##' performs QR decomposition on the data matrices and proceeds from
##' there. However, it is also possible to find the canonical
##' correlations and associated weights using only the covariance
##' matrix. This is especially useful if one wishes to solve the
##' "population" version of CCA for a distribution with population
##' covariance Sigma. This algorithm results in canonical variates
##' that have unit empirical variance.
##'
##' Note: As a convenience to the user, this function will also fix
##' the signs of the resulting coefficients as described in
##' cancor.signfix.
##' @title Canonical Correlation Analysis for Covariance Matrix
##' @param Sigma A PSD matrix given the covariance (sample or
##'     population)
##' @param px The number of variables in the x dataset
##' @return A list containing the following components:
##'
##' cor: correlations
##'
##' xcoef: estimated coefficients for the x variables
##'
##' ycoef: estimated coefficients for the y variables
##' @author Daniel Kessler
##' @export
cancor.cov <- function(Sigma, px) {
  p <- nrow(Sigma)
  sxx <- Sigma[1:px, 1:px]
  syy <- Sigma[(px + 1):p, (px + 1):p]
  syx <- Sigma[(px + 1):p, 1:px]

  sxx.sqrti <- solve(expm::sqrtm(sxx))
  syy.sqrti <- solve(expm::sqrtm(syy))

  SVD <- svd(syy.sqrti %*% syx %*% sxx.sqrti)

  rho <- SVD$d
  xcoef <- sxx.sqrti %*% SVD$v
  ycoef <- syy.sqrti %*% SVD$u

  svd(Sigma)
  fm <- list(cor = rho, xcoef = xcoef, ycoef = ycoef)
  fm <- cancor_signfix(fm)
  return(fm)
}

##' Return the (signed) value of the maximum (in magnitude) element of a vector
##'
##' An example is especially illustrative. absmax(c(1,-3,2)) will
##' yield -3. absmax(c(1,-2,3)) will yield 3.
##' @title Find Maximum Magnitude Element of Vector
##' @param x
##' @return the signed maximum
##' @author Daniel Kessler
absmax <- function(x) {
  x[which.max(abs(x))]
}
