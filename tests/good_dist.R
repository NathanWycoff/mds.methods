#!/usr/bin/Rscript
#  good_dist.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.18.2018

##Test that good_dist does some reasonable things sometimes.
source('../R/forward_methods.R')

# In particular, it should do the same thing for symmetric dissimilarities whether or not we use that fact in computation
n <- 100
p <- 30
X <- matrix(rnorm(n*p), ncol = p)
D1 <- good.dist(X, euclidean.dist, symm = FALSE)
D2 <- good.dist(X, euclidean.dist, symm = TRUE)
print(identical(D1, D2))
