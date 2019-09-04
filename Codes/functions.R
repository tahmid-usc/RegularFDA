## Estimate hypoerparameter


Hyper <- function(dtf) {
  
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
  
  hyp <- optim(par=rep(.01, 5), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par^2)
  
}



# Fit a smooth curve--------------------------------

gpsmooth <- function(y, trainlist) {
  
  ns <- length(y)
  x <- seq(0,1, length.out = trainlist$m)
  k <- covker(y, x, trainlist$theta)
  kxx <- ker1(y, trainlist$theta)
  pred <- k %*% trainlist$mumat %*% trainlist$ySum
  sigmastar <- kxx + trainlist$theta[5] * diag(ns) - trainlist$n * (k %*% trainlist$mumat %*% t(k))
  sigmastar <- forceSymmetric(sigmastar)
  sigmastar <- as.matrix(sigmastar)
  ll <- pred - 1.96 * sqrt(diag(sigmastar))
  ul <- pred + 1.96 * sqrt(diag(sigmastar))
  return(list(pred = pred, sigmastar = sigmastar, ll = ll, ul = ul))
  
}


# Posterior mean and cov matrix for test curve

logprob <- function(fit, testy) {
  
  logprob <- dmvnorm(x = testy, mean = fit$pred, sigma = fitsigmastar, log = T)

  return(logprob)
}


# Extract relevant feature of the data

feature <- function(train) {
  
  n <- dim(train)[1]
  m <- dim(train[,2:19])[2]
  x <- seq(0,1, length.out = m)
  
  #theta <- Hyper(train[,2:19])
  theta <- Hyper.ms(train[,2:19])
  
  k1 <- ker1(x, theta)
  k2 <- ker2(x, theta)
  
  k1inv <- chol2inv(chol(k1))
  D <- chol2inv(chol(k2 + theta[5] * diag(m)))
  C <- chol2inv(chol(k1inv + n * D))
  
  ySum <- as.matrix(apply(train[,2:19], 2, sum), ncol = 1)
  
  mumat <- ( diag(m) - n * D %*% C) %*% D #%*% ySum
  
  return(list(n = n,m = m, theta = theta, mumat = mumat, ySum = ySum))
  
}