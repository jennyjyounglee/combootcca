## Naming conventions used in this file:
## cca_align_*: functions for aligning cca solutions
## cca_ci_*: functions that take data and return CI's (after alignment)
## bt_{prob,algo}_*method*_{inner,fun}: functions for use with batchtools.
##     prob for problems, algo for algorithms.
##     *method* names the approach.
##     fun is for function to be used in prob/algo definition, inner for inrep function




##' Correct sign ambiguity in canonical correlation analysis by requiring that
##' the diagonal of xcoef be non-negative (this is the approach suggested by
##' Anderson 1999). This will also truncate both xcoef and ycoef to have the
##' same number of columns as the length of cor
##'
##' @title Fix cancor signs based on the diagonal
##' @param fm A fitted object returned by cancor (or similar)
##' @param ref Not used
##' @return Object like fm but with possible sign flips
##' @author Dan Kessler
cca_align_posdiag <- function(fm, ref) {
  k <- length(fm$cor)
  fm$xcoef <- fm$xcoef[, 1:k, drop = FALSE]
  fm$ycoef <- fm$ycoef[, 1:k, drop = FALSE]
  signs <- diag(x = sign(diag(fm$xcoef)), nrow = k, ncol = k)
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
cca_align_posmax <- function(fm) {
  K <- min(ncol(fm$xcoef), ncol(fm$ycoef))

  theta <- rbind(fm$xcoef[, 1:K], fm$ycoef[, 1:K])

  maxes <- apply(theta, 2, absmax)

  flip <- diag(sign(maxes))

  fm$xcoef[, 1:K] <- fm$xcoef[, 1:K] %*% flip
  fm$ycoef[, 1:K] <- fm$ycoef[, 1:K] %*% flip

  return(fm)
}

cca_align_greedy_cosx <- function(fm, ref) {
  t_hat <- fm$xcoef
  t_ref <- ref$xcoef

  sim <- cos_sim(t_ref, t_hat)

  ## decompose into signed similarity
  sim_sign <- sign(sim)
  sim_abs <- abs(sim)

  k <- ncol(t_ref)
  map <- rep(0, k)
  signs <- rep(0, k)
  for (i in 1:k) {
    map[i] <- which.max(sim_abs[i, ])
    signs[i] <- sim_sign[i, map[i]]
    sim_abs[, map[i]] <- 0 # zero-out the selected column
  }

  fm$xcoef <- fm$xcoef[, map] %*% diag(signs)
  fm$ycoef <- fm$ycoef[, map] %*% diag(signs)

  return(fm)
}

##' @title Compute cosine similarity between columns of two matrices
##' @param x Matrix of size n by p
##' @param y Matrix of size n by q
##' @return Matrix of size p by q of cosine similarities. The (i, j) entry has
##'   the similarity between x_i and y_j.
##' @author Daniel Kessler
cos_sim <- function(x, y) {
  x <- scale(x, center = FALSE, scale = sqrt(colSums(x^2)))
  y <- scale(y, center = FALSE, scale = sqrt(colSums(y^2)))

  return(t(x) %*% y)
}

##' Compute the canonical correlations between two data matrices such that
##' canonical variates have unit variance.
##'
##' This is a simple wrapper around stats::cancor. stats::cancor uses the Golub
##' and Van Loan algorithm, which yields canonical variates that have unit norm
##' rather than unit variance.
##' As a result, the canonical variates will have (empirical) variance
##' of 1/(N-1), where N is the number of observations. This is
##' undesirable if we wish to perform inference on the weights in
##' xcoef or ycoef, because their scale will vary with N.
##'
##' This function simply calls stats::cancor and then multiplies both xcoef and
##' ycoef by sqrt(N-1) so that the resulting canonical variates will indeed have
##' unit variance.
##'
##' See the documentation for stats::cancor for more details on standard usage.
##'
##' @title Canonical Correlations with Unit Variance
##' @param x one dataset
##' @param y the other dataset
##' @param xcenter whether to center x
##' @param ycenter whether to center y
##' @return
##' @author Daniel Kessler
cancor_scaled <- function(x, y, xcenter = TRUE, ycenter = TRUE) {
  n <- nrow(x)
  fm <- stats::cancor(x, y, xcenter, ycenter)
  fm$xcoef <- sqrt(n - 1) * fm$xcoef
  fm$ycoef <- sqrt(n - 1) * fm$ycoef
  return(fm)
}

##' Obtain confidence intervals for the "directions" of a canonical correlation
##' analysis using asymptotic results from Anderson 1999.
##'
##' Important Note: Theory only valid when p = q.
##' @title Asymptotic confidence intervals for CCA directions
##' @param x Data matrix of size n by p
##' @param y Data matrix of size n by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param align A function to perform post-processing on the estimated
##'   coefficients to render the solution well-identified. By default, this uses
##'   cancor_signfix_diag, which ensures that the diagonal of xcoef is
##'   non-negative. Should also take "ref" or "..." as an argument (but doesn't
##'   have to use it).
##' @param ref A reference solution to align against
##' @return List with two objects: xcoef_ci and ycoef_ci.
##' @author Dan Kessler
##' @export
cca_ci_asymptotic <- function(x, y, level = .95,
                              align = cca_align_posdiag, ref) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  K <- min(p, q)

  fm <- cancor_scaled(x, y)
  fm <- align(fm, ref)
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

##' Bootstrap absolute quantiles to generate CCA confidence intervals
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
##' @return List with two objects: xcoef_ci and ycoef_ci.
##' @author Daniel Kessler
##' @export
cca_ci_absboot <- function(x, y, level = .95, align = cca_align_posdiag, ref,
                             nboots = 1e3, parametric = FALSE, progress = 0) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  fm <- cancor_scaled(x, y)
  fm <- align(fm, ref)
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
      fm_boot <- cancor_scaled(x[idx, ], y[idx, ])
    } else {
      if (i == 1) { # only compute cholesky once
        sigma_hat <- cov(cbind(x, y))
        l_hat <- chol(sigma_hat)
      }
      newdata <- t(l_hat %*% matrix(rnorm(n * (p + q)), p + q, n))
      x_boot <- newdata[, 1:p]
      y_boot <- newdata[, (p + 1):(p + q)]
      fm_boot <- cancor_scaled(x_boot, y_boot)
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
cca_ci_regression <- function(x, y, level = .95, align = cca_align_posdiag, ref,
                              train_ratio = 0.5) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  n_train <- floor(n * train_ratio)
  train_ind <- sample.int(n, size = n_train, replace = FALSE)

  x1 <- x[+train_ind, ]
  x2 <- x[-train_ind, ]
  y1 <- y[+train_ind, ]
  y2 <- y[-train_ind, ]

  fm1 <- cancor_scaled(x1, y1)
  fm1 <- align(fm1, ref)
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

##' Use the boot package to get bootstrapped confidence intervals for CCA
##'
##' Details go here
##' @title Bootstrap-based confidence intervals for CCA
##' @param x Data matrix of size n by p
##' @param y Data matrix of size n by q
##' @param level Level for confidence intervals, should be in (0, 1)
##' @param align Function for alignment
##' @param ref Passed through to align function
##' @param nboots How many bootstrap replicates
##' @return List of several types of CIs
##' @author Dan Kessler
cca_ci_boot <- function(x, y, level=0.90, align = cca_align_posdiag,
                        ref, nboots = 1e2) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  k <- min(p, q, n)

  data_flat <- cbind(x, y)

  boot_stat <- function(data, idx, p, align, ref) {
    bdata <- data[idx, ]
    theta <- cancor_vec(bdata, p, align, ref)
    return(theta)
  }

  boot_out <- boot::boot(data_flat,
    R = nboots,
    statistic = boot_stat,
    p = p,
    align = align,
    ref = ref
  )

  ## preallocate for results
  numel <- (p + q) * k
  ci_norm_flat <- array(NA, c(numel, 2))
  ci_basic_flat <- array(NA, c(numel, 2))
  ci_perc_flat <- array(NA, c(numel, 2))
  ci_bca_flat <- array(NA, c(numel, 2))


  ## loop over theta to get confidence intervals
  for (i in seq_len(numel)) {
    bootci <- boot::boot.ci(
      boot.out = boot_out,
      conf = level,
      type = c("norm", "basic", "perc", "bca"),
      index = i
    )

    ci_norm_flat[i, ] <- bootci$normal[1, 2:3]
    ci_basic_flat[i, ] <- bootci$basic[1, 4:5]
    ci_perc_flat[i, ] <- bootci$percent[1, 4:5]
    ci_bca_flat[i, ] <- bootci$bca[1, 4:5]
  }

  ci_glue <- function(ci_flat) {
    alpha <- 1 - level
    ci_levels <- paste0(c(100 * alpha / 2, 100 * (1 - alpha / 2)), "%")
    adimnames <- list(
      coordinate = NULL,
      component = 1:k,
      ci_levels
    )

    ci_lower <- vec2fm(ci_flat[, 1], p, q)
    ci_upper <- vec2fm(ci_flat[, 2], p, q)

    xcoef_ci <- abind::abind(ci_lower$xcoef, ci_upper$xcoef, along = 3)
    dimnames(xcoef_ci) <- adimnames

    ycoef_ci <- abind::abind(ci_lower$ycoef, ci_upper$ycoef, along = 3)
    dimnames(ycoef_ci) <- adimnames

    fm <- list(xcoef_ci = xcoef_ci, ycoef_ci = ycoef_ci)
    return(fm)
  }

  res <- list(
    ci_norm = ci_glue(ci_norm_flat),
    ci_basic = ci_glue(ci_basic_flat),
    ci_perc = ci_glue(ci_perc_flat),
    ci_bca = ci_glue(ci_bca_flat)
  )

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

  pop.fm <- cancor_cov(Sigma, px)

  boot.cis <- bootstrapcca(X, Y)
}

##' @title Conduct asymptotic coverage experiment
##' @param outreps Each "outer replication" draws a new value of sigma. This can
##'   either be an integer (e.g., 10L), in which case a random covariance matrix
##'   will be drawn for each out replication, or instead it can be a list whose
##'   length determines the number of replications, and whose values are
##'   covariance matrices to use.
##' @param inreps Each "inner replication" is a repetition with new data for the
##'   same value of sigma
##' @param p The dimension of X
##' @param q The dimension of Y
##' @param n How many datapoints to draw in each sample
##' @param nboots How many bootstrap samples to draw
##' @param sigma Optional. A list of length outreps containing covariance
##'   matrices of size p + q. This matrices will be used for each of the outreps
##'   replications. If not provided, a covariance matrix will be generated for
##'   you for each of the outreps.
##' @return A list of results, containing: (1) truth, a 4D array of true
##'   parameter values, (2) cover, a 5D array of logicals indicating which CIs
##'   cover the truth, (3) length, a 5D array of CI lengths, (4) cis, a 6D array
##'   of the confidence intervals, and (5) sigma, a list of length outreps the
##'   holds the generative Sigma. For the arrays, the dimensions (as applicable)
##'   index (i) coordinates (p+q), (ii) components (K), (iii) inreps, (iv)
##'   outreps, (v) method, (vi) lower/upper confidence bounds.
##' @export
##' @author Dan Kessler
coverage_experiment <- function(outreps = 1L, inreps = 1L, p = 2, q = 2,
                                n = 1000L, nboots = 1000L) {
  K <- min(p, q)
  coord_names <- c(paste0("X_", 1:p), paste0("Y_", 1:q))

  methods <- c(
    "asymptotic",
    "split-regression",
    "bootstrap-abs"
  )

  M <- length(methods) # number of methods

  if (is.integer(outreps)) {
    sigma <- list()
    for (i in 1:outreps) {
      sigma[[i]] <- gen_sigma(p, q)
    }
  } else if (is.list(outreps)) {
    sigma <- outreps
    outreps <- length(sigma)
  }

  cis <- array(NA, c((p + q), K, inreps, outreps, M, 2),
    dimnames = list(
      coordinate = coord_names,
      component = 1:K,
      inreps = 1:inreps,
      outreps = 1:outreps,
      method = methods,
      quantity = c("lower", "upper")
    )
  )

  cover <- array(NA, c((p + q), K, inreps, outreps, M),
    dimnames = list(
      coordinate = coord_names,
      component = 1:K,
      inreps = 1:inreps,
      outreps = 1:outreps,
      method = methods)
  )

  length <- cover

  truth <- array(NA, c((p + q), K, inreps, outreps),
    dimnames = list(
      coordinate = coord_names,
      component = 1:K,
      inreps = 1:inreps,
      outreps = 1:outreps)
  )



  for (i in 1:outreps) {
    fm_true <- cancor_cov(sigma[[i]], px = p)
    for (j in 1:inreps) {
      dat <- gen_data(sigma[[i]], p, q, n)
      truth[1:p, , j, i] <- fm_true$xcoef
      truth[(p + 1):(p + q), , j, i] <- fm_true$ycoef

      ## asymptotic
      ci_estimates <- cca_ci_asymptotic(dat$x, dat$y)
      cis[1:p, , j, i, 1, 1] <- ci_estimates$xcoef[, , 1]
      cis[1:p, , j, i, 1, 2] <- ci_estimates$xcoef[, , 2]
      cis[(p + 1):(p + q), , j, i, 1, 1] <- ci_estimates$ycoef[, , 1]
      cis[(p + 1):(p + q), , j, i, 1, 2] <- ci_estimates$ycoef[, , 2]

      ## regression
      ci_estimates <- cca_ci_regression(dat$x, dat$y)
      cis[1:p, , j, i, 2, 1] <- ci_estimates$xcoef[, , 1]
      cis[1:p, , j, i, 2, 2] <- ci_estimates$xcoef[, , 2]
      cis[(p + 1):(p + q), , j, i, 2, 1] <- ci_estimates$ycoef[, , 1]
      cis[(p + 1):(p + q), , j, i, 2, 2] <- ci_estimates$ycoef[, , 2]

      ## bootstrap-abs
      ci_estimates <- cca_ci_bootstrap(dat$x, dat$y, nboots = 10)
      cis[1:p, , j, i, 3, 1] <- ci_estimates$xcoef[, , 1]
      cis[1:p, , j, i, 3, 2] <- ci_estimates$xcoef[, , 2]
      cis[(p + 1):(p + q), , j, i, 3, 1] <- ci_estimates$ycoef[, , 1]
      cis[(p + 1):(p + q), , j, i, 3, 2] <- ci_estimates$ycoef[, , 2]
    }
  }

  for (m in 1:M) {
    cover[, , , , m] <- cis[, , , , m, 1] <= drop(truth) &
      drop(truth) <= cis[, , , , m, 2]

    length[, , , , m] <- cis[, , , , m, 2] - cis[, , , , m, 1]
  }

  return(list(cis = cis, truth = truth, cover = cover,
              length = length, sigma = sigma))
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
cancor_cov <- function(Sigma, px, align = cca_align_posdiag) {
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

## vectorized version of cancor for use with boot
cancor_vec <- function(data, p, align, ref) {
  q <- ncol(data) - p
  x <- data[, 1:p]
  y <- data[, (p + 1):(p + q)]
  fm <- cancor_scaled(x, y)
  fm <- align(fm, ref)
  theta <- fm2vec(fm)
  return(theta)
}

## batchtools functions below here

##' This is a convenience function for use with batchtools
##'
##' The details go here
##' @title Batchtools convenience function for "problem"
##' @param job Constructed by batchtools
##' @param data Optionally, pass a list with named elements `x` and `y`, where
##'   each is an n by p and n by q matrix, respectively
##' @param sigma Optional: a square covariance matrix of dimension (p + q)
##' @param p The dimension of the random variable X
##' @param q The dimension of the random variable Y
##' @param n The number of observations
##' @param inreps How many replicates to draw
##' @return An instance
##' @author Daniel Kessler
##' @export
bt_problem_std_fun <- function(job = NULL, data = NULL, sigma = NULL, p = NULL,
                               q = NULL, n = NULL, inreps = 1L) {
  prob_fun_inner <- function(data = NULL) {
    res <- list()
    if (!is.null(data)) res$data <- data
    if (is.null(data)) res$data <- ccasleuth:::gen_data(sigma, p, q, n)
    res$fm_hat <- cancor_scaled(res$data$x, res$data$y)
    return(res)
  }

  instance <- list()
  instance$inreps <- list()

  if (!is.null(data)) instance$inreps[1] <- prob_fun_inner(data = data)

  if (is.null(data)) {
    if (is.null(sigma)) sigma <- ccasleuth:::gen_sigma(p, q)
    instance$inreps <- replicate(inreps, prob_fun_inner(), simplify = FALSE)
  }

  instance$sigma <- sigma
  if (!is.null(sigma)) instance$fm_true <- cancor_cov(sigma, p)

  return(instance)
}

##' @title Generic inrep function for batchtools algorithms
##' @param inrep The inrep object
##' @param ci_func A function that follows the cca_ci_* API
##' @param ... Arguments passed on to ci_func
##' @return CIs as returned by ci_func
##' @author Dan Kessler
bt_algo_inrep <- function(inrep, ci_func, ...) {
  res <- ci_func(
    x = inrep$data$x,
    y = inrep$data$y,
    ref = inrep$fm_hat,
    ...
  )
  return(res)
}

##' @title Batchtools Function for Asymptotic Confidence Intervals
##' @param job required by batchtools
##' @param data required by batchtools
##' @param instance required by batchtools
##' @param ... Arguments passed on to ci_func
##' @return CI object
##' @export
##' @author Dan Kessler
bt_algo_asymptotic <- function(job, data, instance, ...) {
  res <- lapply(instance$inreps, bt_algo_inrep,
    ci_func = cca_ci_asymptotic,
    align = cca_align_posdiag,
    ...
  )
  return(res)
}

##' @title Batchtools Function for Bootstrapped Abs
##' @param job required by batchtools
##' @param data required by batchtools
##' @param instance required by batchtools
##' @param ... Arguments passed on to ci_func
##' @return CI object
##' @export
##' @author Dan Kessler
bt_algo_absboot <- function(job, data, instance, ...) {
  res <- lapply(instance$inreps, bt_algo_inrep,
    ci_func = cca_ci_absboot,
    align = cca_align_posdiag,
    ...
  )
  return(res)
}

##' @title Batchtools Function for Regression-Based CIs
##' @param job required by batchtools
##' @param data required by batchtools
##' @param instance required by batchtools
##' @param ... Arguments passed on to ci_func
##' @return CI object
##' @export
##' @author Dan Kessler
bt_algo_regression <- function(job, data, instance, ...) {
  res <- lapply(instance$inreps, bt_algo_inrep,
    ci_func = cca_ci_regression,
    align = cca_align_posdiag,
    ...
  )
  return(res)
}

##' @title Batchtools Function for Bootstrapped Confidence Intervals
##' @param job required by batchtools
##' @param data required by batchtools
##' @param instance required by batchtools
##' @param ... Arguments passed on to ci_func
##' @return CI object
##' @export
##' @author Dan Kessler
bt_algo_boot <- function(job, data, instance, ...) {
  res <- lapply(instance$inreps, bt_algo_inrep,
    ci_func = cca_ci_boot,
    align = cca_align_posdiag,
    ...
  )
  return(res)
}


fm2mat <- function(fm) rbind(fm$xcoef, fm$ycoef)

fm2vec <- function(fm) mat2vec(fm2mat(fm))

mat2fm <- function(mat, p) {
  fm <- list()
  fm$xcoef <- mat[1:p, ]
  fm$ycoef <- mat[(p + 1):nrow(mat), ]
  return(fm)
}

mat2vec <- function(mat) c(mat)

vec2fm <- function(vec, p, q) mat2fm(vec2mat(vec, p, q), p)

vec2mat <- function(vec, p, q) matrix(vec, nrow = p + q)
