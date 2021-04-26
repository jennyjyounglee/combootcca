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
bootstrapcca <- function(X, Y, level = .05, nboots = 1e3){
    N <- nrow(X)
    px <- ncol(X)
    py <- ncol(Y)
    K <- min(px, py)

    fm <- cancor(X, Y)
    fm <- cancor.scalefix(cancor.signfix(fm))
    rho.hat <- fm$cor
    beta.hat <- fm$ycoef
    alpha.hat <- fm$xcoef

    rho.boot <- array(dim=c(length(rho.hat),nboots))
    alpha.boot <- array(dim=c(dim(alpha.hat),nboots))
    beta.boot <- array(dim=c(dim(beta.hat),nboots))

    for (i in 1:nboots){
        idx <- sample(N, replace=TRUE)
        fm.boot <- cancor(X[idx,], Y[idx,])
        fm.boot <- cancor.scalefix(cancor.signfix(fm))
        rho.boot[,i] <- fm.boot$cor
        alpha.boot[,,i] <- fm.boot$xcoef
        beta.boot[,,i] <- fm.boot$ycoef
    }

    rho.dist  <- abs(sweep(rho.boot, c(1), rho.hat))
    alpha.dist <- abs(sweep(alpha.boot, c(1, 2), alpha.hat))
    beta.dist <- abs(sweep(beta.boot, c(1, 2), beta.hat))

    rho.t <- apply(rho.dist, 1, quantile, probs = level)
    alpha.t <- apply(alpha.dist, c(1,2), quantile, probs = level)
    beta.t <- apply(beta.dist, c(1,2), quantile, probs = level)

    rho.ci <- array(c(rho.hat - rho.t, rho.hat, rho.hat + rho.t),
                    dim = c(K, 3),
                    dimnames = list(Component = 1:K, Interval = c("Lower", "Estimate", "Upper")))

    alpha.ci <- array(c(
        alpha.hat - alpha.t,
        alpha.hat,
        alpha.hat + alpha.t),
        dim = c(dim(alpha.hat), 3),
        dimnames = list(Coordinate = 1:px, Component = 1:px, Interval = c("Lower", "Estimate", "Upper")))
    alpha.ci <- alpha.ci[,1:K,] # only worry about K components

    beta.ci <- array(c(
        beta.hat - beta.t,
        beta.hat,
        beta.hat + beta.t),
        dim = c(dim(beta.hat), 3),
        dimnames = list(Coordinate = 1:py, Component = 1:py, Interval = c("Lower", "Estimate", "Upper")))

    beta.ci <- beta.ci[,1:K,] # only worry about K components

    CIs = list(rho = rho.ci, alpha = alpha.ci, beta = beta.ci)
    return(CIs)
    
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
cancor.signfix <- function(fm){
    K <- min(ncol(fm$xcoef), ncol(fm$ycoef))

    theta <- rbind(fm$xcoef[,1:K], fm$ycoef[,1:K])
    
    maxes <- apply(theta, 2, absmax)

    flip <- diag(sign(maxes))

    fm$xcoef[,1:K] <- fm$xcoef[,1:K] %*% flip
    fm$ycoef[,1:K] <- fm$ycoef[,1:K] %*% flip

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
cancor.scalefix <- function(fm, N){
    fm$xcoef <- sqrt(N-1) * fm$xcoef
    fm$ycoef <- sqrt(N-1) * fm$ycoef
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
cancor.cov <- function(Sigma, px){
    p <- nrow(Sigma)
    sxx <- Sigma[1:px,1:px]
    syy <- Sigma[(px+1):p,(px+1):p]
    syx <- Sigma[(px+1):p,1:px]

    sxx.sqrti <- solve(expm::sqrtm(sxx))
    syy.sqrti <- solve(expm::sqrtm(syy))

    SVD <- svd(syy.sqrti %*% syx %*% sxx.sqrti)

    rho <- SVD$d
    xcoef <- sxx.sqrti %*% SVD$v
    ycoef <- syy.sqrti %*% SVD$u

    svd(Sigma)
    fm <- list(cor = rho, xcoef = xcoef, ycoef = ycoef)
    fm <- cancor.signfix(fm)
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
absmax <- function(x){
    x[which.max( abs(x) )]
}
