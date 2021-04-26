## code to support simulation studies of CCA

generate.sigma <- function(px, py){
    p <- px + py
    Sigma <- clusterGeneration::genPositiveDefMat(p)$Sigma
    return(Sigma)
}

generate.data <- function(n, Sigma, px){
    p <- nrow(Sigma)
    Z <- MASS::mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
    X <- Z[,1:px]
    Y <- Z[,(px+1):ncol(Sigma)]
    return(list(X=X,Y=Y))
}

sim.coverage.check <- function(n, px, py){
    Sigma <- generate.sigma(px, py)
    dat <- generate.data(n, Sigma, px)
    X <- dat$X
    Y <- dat$Y

    pop.fm <- cancor.cov(Sigma, px)

    boot.cis <- bootstrapcca(X, Y)

}
