##' Bootstrap data and
##'
##' .. content for \details{} ..
##' @title Bootstrap-based confidence intervals for CCA directions
##' @param x Data matrix of size N by p
##' @param y Data matrix of size N by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param nboots Number of bootstrap sample to draw
##' @param parametric If FALSE (default), do bootstrap by sampling with
##'   replacement. If TRUE, perform parametric bootstrap, i.e., draw data from
##'   multivariate normal distribution following sample covariance
##' @param progress If 0 (default), don't report progress. If set to a positive
##'   integer k, report bootstrap progress every k'th run.
##' @return
##' @author Daniel Kessler
##' @export
cca_ci_bootstrap <- function(x, y, level = .05, nboots = 1e3, parametric = FALSE,
                             progress = 0) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  fm <- cancor_scalefix(cancor(x, y), n)
  K <- length(fm$cor)

  rho_hat <- fm$cor
  xcoef_hat <- fm$xcoef[, 1:K]
  ycoef_hat <- fm$ycoef[, 1:K]

  rho_boot <- array(dim = c(length(rho_hat), nboots))
  xcoef_boot <- array(dim = c(dim(xcoef_hat), nboots))
  ycoef_boot <- array(dim = c(dim(ycoef_hat), nboots))

  for (i in 1:nboots) {
    if (progress && !(i %% progress)) {
      progress_message <- sprintf("Running bootstrap %d of %d", i, nboots)
      print(progress_message)
    }
    idx <- sample(n, replace = TRUE)
    if (!parametric) {
      idx <- sample(n, replace = TRUE)
      fm_boot <- cancor_scalefix(cancor(x[idx, ], y[idx, ]), n)
    } else {
      if (i == 1) { # only compute cholesky once
        sigma_hat <- cov(cbind(x, y))
        l_hat <- chol(sigma_hat)
      }
      newdata <- t(l_hat %*% matrix(rnorm(n * (p + q)), p + q, n))
      x_boot <- newdata[, 1:p]
      y_boot <- newdata[, (p + 1):(p + q)]
      fm_boot <- cancor_scalefix(cancor(x_boot, y_boot), n)
    }
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
  xcoef_ci[, , 2] <- xcoef_hat + xcoef_t

  ycoef_ci[, , 1] <- ycoef_hat - ycoef_t
  ycoef_ci[, , 2] <- ycoef_hat + ycoef_t


  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}
