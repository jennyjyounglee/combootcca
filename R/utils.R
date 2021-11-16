## utility functions common to multiple routines


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
