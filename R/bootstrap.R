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

    fm <- cancor.signfix(X, Y)
    rho.hat <- fm$cor
    beta.hat <- fm$ycoef
    alpha.hat <- fm$xcoef

    rho.boot <- array(dim=c(length(rho.hat),nboots))
    alpha.boot <- array(dim=c(dim(alpha.hat),nboots))
    beta.boot <- array(dim=c(dim(beta.hat),nboots))

    for (i in 1:nboots){
        idx <- sample(N, replace=TRUE)
        fm.boot <- cancor.signfix(X[idx,], Y[idx,])
        rho.boot[,i] <- fm.boot$cor
        alpha.boot[,,i] <- fm.boot$xcoef
        beta.boot[,,i] <- fm.boot$ycoef
    }

    rho.dist  <- abs(sweep(rho.boot, c(1), rho.hat))
    alpha.dist <- abs(sweep(alpha.boot, c(1, 2), alpha.hat))
    beta.dist <- abs(sweep(beta.boot, c(1, 2), beta.hat))
    
}

##' Compute canonical correlations between two data matrices, but
##' impose identifiability conditions so the signs are fixed rather
##' than random
##'
##' Canonical Correlation Analysis suffers from sign ambiguity in the
##' estimated coefficients. This is because Corr(U, V) = Corr(-U, -V)
##' and relatedly obtains from the SVD decomposition typically used to
##' solve the CCA problem. This function adopts the convention that
##' each canonical coefficient pair (xcoef, ycoef) should satisfy the
##' condition that its maximal element (in absolute value) has
##' positive sign. You may wish to standardize the columns of X and Y
##' prior to calling this function so that differences in variable
##' scaling do not dominate the identity of the maximal element.
##' @title Sign-Fixed Canonical Correlations
##' @param ... Arguments passed on to cancor
##' @return Same object as returned by cancor after sign-flipping per
##'     the identifiability condition discussed in Details.
##' @author Daniel Kessler
cancor.signfix <- function(...){
    fm <- cancor(...)

    K <- min(ncol(fm$xcoef), ncol(fm$ycoef))

    theta <- rbind(fm$xcoef[,1:K], fm$ycoef[,1:K])
    
    maxes <- apply(theta, 2, absmax)

    flip <- diag(sign(maxes))

    fm$xcoef[,1:K] <- fm$xcoef[,1:K] %*% flip
    fm$ycoef[,1:K] <- fm$ycoef[,1:K] %*% flip

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
