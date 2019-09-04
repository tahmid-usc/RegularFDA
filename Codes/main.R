#
source('code/RBF.R')
source('code/matern52.R')
source('code/laplace.R')
source('code/matern32.R')

source('code/functions.R')


fet <- lapply(train, feature)

y <- seq(0,1, length.out = 100)
mu1 <- gpsmooth(y, fet$G1)
mu2 <- gpsmooth(y, fet$NonG1)


par(mfrow = c(2,1))
plot(x, train$G1[1,2:19], ylim = c(min(train$G1[,2:19]), max(train$G1[,2:19])), col = 'gray30', cex = .6)
apply(train$G1[,2:19], 1, function(y) {points(x, y, col = 'gray30', cex = .6)})
lines(y, mu1$pred, lwd = 2)
lines(y, mu1$ul, lwd = 2, lty = 2)
lines(y, mu1$ll, lwd = 2, lty = 2)
plot(x, train$NonG1[1,2:19], ylim = c(min(train$NonG1[,2:19]), max(train$NonG1[,2:19])), col = 'gray30', cex = .6)
apply(train$NonG1[,2:19], 1, function(y) {points(x, y, col = 'gray30', cex = .6)})
lines(y, mu2$pred, lwd = 2)
lines(y, mu2$ul, lwd = 2, lty = 2)
lines(y, mu2$ll, lwd = 2, lty = 2)


y <- seq(0,1, length.out = 18)
fit1 <- gpsmooth(y, fet$G1)
fit2 <- gpsmooth(y, fet$NonG1)

log.prob <- c()
for (i in 1:dim(mdt)[1]) {
  log.prob1 <- dmvnorm(x = mdt[i, 2:19], mean = fit1$pred, sigma = fit1$sigmastar, log = T) + log(222/612)
  log.prob2 <- dmvnorm(x = mdt[i, 2:19], mean = fit2$pred, sigma = fit2$sigmastar, log = T) + log(390/612)
  log.prob <- rbind(log.prob, cbind(log.prob1, log.prob2))
}

classpred <- ifelse(log.prob[,1] > log.prob[,2], 'G1', 'NonG1')
table(Pred = classpred, True = mdt$phase)


