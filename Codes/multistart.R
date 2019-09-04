#
library(optimx)

Hyper.ms <- function(dtf) {
  
  marlik <- function(theta) {
    
    n <- dim(dtf)[1]
    m <- dim(dtf)[2]
    x <- seq(0,1, length.out = m)
    theta <- theta^2
    
    k1 <- ker1(x, theta)
    k2 <- ker2(x, theta)
    
    k1inv <- chol2inv(chol(k1))
    G <- k2 + theta[5] * diag(m)
    D <- chol2inv(chol(G))
    C <- chol2inv(chol(k1inv + n * D))
    ySum <- as.matrix(apply(dtf, 2, sum), ncol = 1)
    P <- D %*% ySum
    YD <- sum(apply(dtf, 1, function(x) { x <- matrix(x, ncol = 1); t(x) %*% D %*% x}))
    logl <- YD - (t(P) %*% C %*% P) + determinant((k1inv + n * D), logarithm = T)$modulus + determinant(k1, logarithm = T)$modulus + 
      n * determinant(G, logarithm = T)$modulus
    
    return(logl)
  }
  
  parval <- c(.01, .05, .1, .5, 1, 1.5, 2, 3, 5)
  parmat <- matrix(parval, nrow = length(parval), ncol = 5)
  hyp <- multistart(parmat, fn = marlik, method='Nelder-Mead', control = list(maxit = 10000))
  hyp <- as.numeric(hyp[which(hyp$value == min(hyp$value)), 1:5])^2
  return(hyp)
  
}