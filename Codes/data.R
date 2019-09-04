rm(list=ls())
library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(readxl)

yeast <- read.delim("data/yeast.txt")
gene <- read_excel("data/gene.xlsx")
names(gene)[1] <- "X"
mdt <- merge(yeast, gene, by = "X")
phase <- ifelse(mdt$Peak == 'G1', 'G1', 'NonG1')
mdt <- cbind(mdt, phase)
mdt <- mdt[,-(2:7)]
mdt <- mdt[, -(20:78)]
mdt <- na.omit(mdt)

sel <- sample(1:dim(mdt)[1], 300)

##mdt <- mdt[sel,]

train <- split(mdt, mdt$phase, drop = T)

x <- seq(0,1, length.out = 18)

par(mfrow = c(2,1))
plot(x, train$G1[1,2:19], type = 'b', ylim = c(min(train$G1[,2:19]), max(train$G1[,2:19])), col = 'gray30')
apply(train$G1[,2:19], 1, function(y) {lines(x, y, type = 'b', col = 'gray30')})
plot(x, train$G1[1,2:19], type = 'b', ylim = c(min(train$G1[,2:19]), max(train$NonG1[,2:19])), col = 'gray30')
apply(train$NonG1[,2:19], 1, function(y) {lines(x, y, type = 'b', col = 'gray30')})

