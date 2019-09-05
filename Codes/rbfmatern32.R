## RBF kernel matern 32

ker1 <- function(x, theta, jitter = .00000001) {
  m <- length(x)
  k1ker <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(k1ker, x = x)
  k1 <- k1 + jitter * diag(m)
  return(k1)
}


ker2 <- function(x, theta) {
  m <- length(x)
  r <- as.matrix(dist(x))
  k2 <- theta[4] * (1 + (sqrt(3) * r / theta[3])) * exp(- sqrt(3) * r /theta[3])
  return(k2)
}

covker <- function(x, y, theta) {
  rbf <- rbfdot(sigma = 1/theta[1])
  covker <- theta[2] * kernelMatrix(rbf, x = x, y = y)
  return(covker)
}

