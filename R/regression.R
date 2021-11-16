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


  xcoef_ci <- array(NA, c(p, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
  )

  ycoef_ci <- array(NA, c(q, K, 2),
    dimnames = list(NULL, NULL, paste(c(100 * (1 - level), 100 * level), "%"))
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
