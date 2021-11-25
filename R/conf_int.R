##' Obtain confidence intervals for the "directions" of a canonical correlation
##' analysis using asymptotic results from Anderson 1999.
##'
##' @title Asymptotic confidence intervals for CCA directions
##' @param x Data matrix of size n by p
##' @param y Data matrix of size n by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param align A function to perform post-processing on the estimated
##'   coefficients to render the solution well-identified. By default, this uses
##'   cancor_signfix_diag, which ensures that the diagonal of xcoef is
##'   non-negative.
##' @return List with two objects: xcoef_ci and ycoef_ci.
##' @author Dan Kessler
##' @export
cca_ci_asymptotic <- function(x, y, level = .95, align = cancor_signfix_diag) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  K <- min(p, q)

  fm <- cancor(x, y) # fit the CCA model
  fm <- align(cancor_scalefix(fm, n))
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

  alpha <- 1 - level
  ci_levels <- paste(c(100 * alpha / 2, 100 * (1 - alpha / 2)), "%")

  xcoef_ci <- array(NA, c(p, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  ycoef_ci <- array(NA, c(q, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  xcoef_ci[, , 1] <- fm$xcoef - sqrt(xvar / n) * zcrit
  xcoef_ci[, , 2] <- fm$xcoef + sqrt(xvar / n) * zcrit


  ycoef_ci[, , 1] <- fm$ycoef - sqrt(yvar / n) * zcrit
  ycoef_ci[, , 2] <- fm$ycoef + sqrt(yvar / n) * zcrit


  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}

##' Bootstrap data to generate CCA confidence intervals
##'
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

  alpha <- 1 - level
  ci_levels <- paste(c(100 * alpha / 2, 100 * (1 - alpha / 2)), "%")

  xcoef_ci <- array(NA, c(p, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  ycoef_ci <- array(NA, c(q, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  xcoef_ci[, , 1] <- xcoef_hat - xcoef_t
  xcoef_ci[, , 2] <- xcoef_hat + xcoef_t

  ycoef_ci[, , 1] <- ycoef_hat - ycoef_t
  ycoef_ci[, , 2] <- ycoef_hat + ycoef_t


  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}

##' Obtain confidence intervals for canonical correlation analysis directions by
##' taking a regression-approach
##'
##' The idea of the regression approach is to recast the CCA problem as a
##' regression problem. First, select a proportion of the data equal to
##' \code{train_ratio} of the data (rounding down as necessary). Then, estimate
##' the usual CCA parameters on this partition: these include \code{cor} (the
##' canonical correlation coefficients), \code{xcoef}, (the canonical directions
##' associated with x), and \code{ycoef} (the canonical directions associated
##' with y).
##'
##' Use the canonical directions estimated above to construct new canonical
##' variables in the held-out partition (by taking the inner product between
##' each observation and a given canonical direction). Then, fit many regression
##' models where we predict one canonical variable using the other dataset. For
##' each, apply a transformation such that the predicted values have unit
##' variance (to satisfy the CCA constraint). Finally, use conventional
##' regression-based inference to provide confidence intervals for the
##' coefficients.
##' @title Regression-Based Confidence Intervals for CCA Directions
##' @param x Data matrix of size n by p
##' @param y Data matrix of size n by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param train_ratio What proportion of the data should be used for estimating
##'   CCA directions
##' @return List with two objects: xcoef_ci and ycoef_ci.
##' @author Dan Kessler
##' @export
cca_ci_regression <- function(x, y, level = .95, train_ratio = 0.5) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  n_train <- floor(n * train_ratio)
  train_ind <- sample.int(n, size = n_train, replace = FALSE)

  x1 <- x[+train_ind, ]
  x2 <- x[-train_ind, ]
  y1 <- y[+train_ind, ]
  y2 <- y[-train_ind, ]

  fm1 <- cancor_scalefix(cancor(x1, y1), n_train)
  K <- length(fm1$cor)

  x2scores <- x2 %*% fm1$xcoef[, 1:K]
  y2scores <- y2 %*% fm1$ycoef[, 1:K]

  alpha <- 1 - level
  ci_levels <- paste(c(100 * alpha / 2, 100 * (1 - alpha / 2)), "%")

  xcoef_ci <- array(NA, c(p, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  ycoef_ci <- array(NA, c(q, K, 2),
                    dimnames = list(NULL, NULL, ci_levels)
  )

  for (k in 1:K) {
    fm2_xpred <- lm_scalefix(lm(y2scores[, k] ~ x2))
    fm2_ypred <- lm_scalefix(lm(x2scores[, k] ~ y2))


    xcoef_ci[, k, 1] <- confint(fm2_xpred)[-1, 1]
    xcoef_ci[, k, 2] <- confint(fm2_xpred)[-1, 2]

    ycoef_ci[, k, 1] <- confint(fm2_ypred)[-1, 1]
    ycoef_ci[, k, 2] <- confint(fm2_ypred)[-1, 2]
  }

  res <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
  return(res)
}

##' @title Refit a linear model such that predicted values have unit variance
##' @param fm An object fitted with lm
##' @return An object of type lm refit such that predicted values have unit
##'   variance
##' @author Dan Kessler
lm_scalefix <- function(fm) {
  y <- fm$model[[1]] # extract response
  x <- fm$model[[2]] # extract predictors
  y_tilde <- y / sd(predict(fm))
  fm_tilde <- lm(y_tilde ~ x)
  return(fm_tilde)
}
