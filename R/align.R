## functions related to aligning CCA solutions

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
