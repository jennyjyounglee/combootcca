##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Asymptotic Confidence Intervals for CCA Directions
##' @param x Data matrix of size n by p
##' @param y Data matrix of size n by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @return List with two objects: xcoef_ci and ycoef_ci.
##' @author Dan Kessler
cca_ci_asymptotic <- function(x, y, level = .95) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  K <- min(p, q)

  fm <- cancor(x, y) # fit the CCA model
  fm <- cancor_scalefix(fm, n)
  fm$xcoef <- fm$xcoef[, 1:K]
  fm$ycoef <- fm$ycoef[, 1:K]

  xvar <- matrix(NA, nrow = nrow(fm$xcoef), ncol = ncol(fm$xcoef))
  yvar <- matrix(NA, nrow = nrow(fm$ycoef), ncol = ncol(fm$ycoef))

  for (j in 1:K) {
    for (i in 1:p) {
      runsum <- 0
      for (k in setdiff(1:K, j)) {
        runsum <- runsum +
          with(
            fm,
            (cor[k]^2 + cor[j]^2 - 2 * cor[k]^2 * cor[j]^2) /
              (cor[j]^2 - cor[k]^2)^2 *
              xcoef[i, k]^2
          )
      }
      xvar[i, j] <- 1 / 2 * fm$xcoef[i, j]^2 + (1 - fm$cor[j]^2) * runsum
    }
  }

  for (j in 1:K) {
    for (i in 1:q) {
      runsum <- 0
      for (k in setdiff(1:K, j)) {
        runsum <- runsum +
          with(
            fm,
            (cor[k]^2 + cor[j]^2 - 2 * cor[k]^2 * cor[j]^2) /
              (cor[j]^2 - cor[k]^2)^2 *
              ycoef[i, k]^2
          )
      }
      yvar[i, j] <- 1 / 2 * fm$ycoef[i, j]^2 + (1 - fm$cor[j]^2) * runsum
    }
  }

    zcrit <- - qnorm((1 - level) / 2)

  xcoef_ci <- array(NA, c(p, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
  )

  ycoef_ci <- array(NA, c(q, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
  )

  xcoef_ci[, , 1] <- fm$xcoef - sqrt(xvar / n) * zcrit
  xcoef_ci[, , 2] <- fm$xcoef + sqrt(xvar / n) * zcrit


  ycoef_ci[, , 1] <- fm$ycoef - sqrt(yvar / n) * zcrit
  ycoef_ci[, , 2] <- fm$ycoef + sqrt(yvar / n) * zcrit


  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}
