##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title bootstrap CCA
##' @param X Data matrix of size N by p
##' @param Y Data matrix of size N by q
##' @param nboots Number of bootstrap sample to draw
##' @return 
##' @author Daniel Kessler
bootstrapcca <- function(X, Y, nboots){
    N <- nrow(X)
    fm <- cancor(X, Y)
    rho.hat <- fm$cor
    beta.hat <- fm$ycoef
    alpha.hat <- fm$xcoef

    rho.boot <- array(dim=c(length(rho.hat),nboots))
    beta.boot <- array(dim=c(dim(beta.hat),nboots))
    alpha.boot <- array(dim=c(dim(alpha.hat),nboots))

    for (i in 1:nboots){
        idx <- sample(N, replace=TRUE)
        fm.boot <- cancor(X[idx,], Y[idx,])
        rho.boot[,i] <- fm.boot$cor
        beta.boot[,,i] <- fm.boot$ycoef
        alpha.boot[,,i] <- fm.boot$xcoef
    }
    
}
