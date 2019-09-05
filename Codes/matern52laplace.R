# #Matern laplace

ker1 <- function(x, theta, jitter = .00000001) {
  m <- length(x)
  r <- as.matrix(dist(x))
  k1 <- theta[2] * (1 + (sqrt(5) * r / theta[1]) + ((5 * r^2) / (3 * theta[1]^2))) * exp(- sqrt(5) * r / theta[1])
  k1 <- k1 + jitter * diag(m)
  return(k1)
}


ker2 <- function(x, theta) {
  k2ker <- laplacedot(sigma = 1/theta[3]) 
  k2 <- theta[4] * kernelMatrix(k2ker, x = x) 
  return(k2)
}


covker <- function(x, y, theta) {
  r <- abs(outer(x, y, "-"))
  covker <- theta[2] * (1 + (sqrt(5) * r / theta[1]) + ((5 * r^2) / (3 * theta[1]^2))) * exp(- sqrt(5) * r / theta[1])
  return(covker)
}


