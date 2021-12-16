## functions related to aligning CCA solutions

##' Correct sign ambiguity in canonical correlation analysis by requiring that
##' the diagonal of xcoef be non-negative (this is the approach suggested by
##' Anderson 1999). This will also truncate both xcoef and ycoef to have the
##' same number of columns as the length of cor
##'
##' @title Fix cancor signs based on the diagonal
##' @param fm A fitted object returned by cancor (or similar)
##' @return Object like fm but with possible sign flips
##' @author Dan Kessler
cancor_signfix_diag <- function(fm) {
  k <- length(fm$cor)
  fm$xcoef <- fm$xcoef[, 1:k]
  fm$ycoef <- fm$ycoef[, 1:k]
  signs <- diag(sign(diag(fm$xcoef)))
  fm$xcoef <- fm$xcoef %*% signs
  fm$ycoef <- fm$ycoef %*% signs
  return(fm)
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
##' @title Fix cancor signs based on the max entry in each column
##' @param fm A fitted object returned by cancor
##' @return Same object as returned by cancor after sign-flipping per
##'     the identifiability condition discussed in Details.
##' @author Daniel Kessler
cancor_signfix_max <- function(fm) {
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
cca_ci_bootstrap <- function(x, y, level = .95, nboots = 1e3, parametric = FALSE,
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
## code to support simulation studies of CCA

sim.coverage.check <- function(n, px, py) {
  Sigma <- generate.sigma(px, py)
  dat <- generate.data(n, Sigma, px)
  X <- dat$X
  Y <- dat$Y

  pop.fm <- cancor.cov(Sigma, px)

  boot.cis <- bootstrapcca(X, Y)
}

##' @title Conduct asymptotic coverage experiment
##' @param outreps Each "outer replication" draws a new value of sigma
##' @param inreps Each "inner replication" is a repetition with new data for the
##'   same value of sigma
##' @param p The dimension of X
##' @param q The dimension of Y
##' @param n How many datapoints to draw in each sample
##' @param ci_method A list of methods for constructing CCA confidence
##'   intervals. Needs an API-like cca_ci_asymptotic (the default)
##' @param sigma Optional. A list of length outreps containing covariance
##'   matrices of size p + q. This matrices will be used for each of the outreps
##'   replications. If not provided, a covariance matrix will be generated for
##'   you for each of the outreps.
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
                                ci_method = cca_ci_asymptotic,
                                sigma, ...) {
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
##' @title Canonical Correlation Analysis for Covariance Matrix
##' @param Sigma A PSD matrix given the covariance (sample or population)
##' @param px The number of variables in the x dataset
##' @param align A function to perform post-processing on the estimated
##'   coefficients to render the solution well-identified. By default, this uses
##'   cancor_signfix_diag, which ensures that the diagonal of xcoef is
##'   non-negative.
##' @return A list containing the following components:
##'
##' cor: correlations
##'
##' xcoef: estimated coefficients for the x variables
##'
##' ycoef: estimated coefficients for the y variables
##' @author Daniel Kessler
##' @export
cancor.cov <- function(Sigma, px, align = cancor_signfix_diag) {
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
  fm <- align(fm)
  return(fm)
}

##' Return the (signed) value of the maximum (in magnitude) element of a vector
##'
##' An example is especially illustrative. absmax(c(1,-3,2)) will
##' yield -3. absmax(c(1,-2,3)) will yield 3.
##' @title Find Maximum Magnitude Element of Vector
##' @param x A vector
##' @return the signed maximum
##' @author Daniel Kessler
absmax <- function(x) {
  x[which.max(abs(x))]
}

## vectorized version of cancor
cancor_vec <- function(data, p, align = cancor_signfix_diag) {
  n <- nrow(data)
  q <- ncol(data) - p
  x <- data[, 1:p]
  y <- data[, (p + 1):(p + q)]
  fm <- align(cancor_scalefix(cancor(x, y), n))
  theta <- c(fm$xcoef, fm$ycoef)
  return(theta)
}
