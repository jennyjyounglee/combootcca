regressioncca <- function(X, Y){
    N <- nrow(X)
    n.train <- floor(N/2)
    train.ind <- sample.int(N, size = N.train, replace = FALSE)

    X1 <- X[+train.ind,]
    X2 <- X[-train.ind,]
    Y1 <- Y[+train.ind,]
    Y2 <- Y[-train.ind,]

    
}
