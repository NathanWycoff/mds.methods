#!/usr/bin/Rscript
#  fast_dist.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.18.2018

euclidean.dist <- function(x, y) norm(x-y, '2')

## How can we speed up distance computation? Try some different ways

## the existing way
standard <- function(x, dist.func, weights = null) {
    if (!is.null(weights)) {
        weights <- sqrt(weights)
        x <- x %*% diag(weights)
    }
    n <- dim(x)[1]
    dists <- matrix(0, ncol = n, nrow = n)
    for (i in 1:n) {
        for (j in 1:n) {
            dists[i,j] <- dist.func(x[i,], x[j,])
        }
    }
    return(dists)
}

## With symmetry
symm <- function(X, dist.func, weights = NULL, symm = FALSE) {
    if (!is.null(weights)) {
        weights <- sqrt(weights)
        X <- X %*% diag(weights)
    }
    n <- dim(X)[1]
    dists <- matrix(0, ncol = n, nrow = n)

    #Calculate the actual distances
    if (symm) {
        for (i in 1:n) {
            for (j in 1:i) {
                dists[i,j] <- dists[j,i] <- dist.func(X[i,], X[j,])
            }
        } 
    } else {
        for (i in 1:n) {
            for (j in 1:n) {
                dists[i,j] <- dist.func(X[i,], X[j,])
            }
        }
    }
    return(dists)
}

##Testing
n <- 300
p <- 30
X <- matrix(rnorm(n*p), ncol = p)

system.time(Ds <- standard(X, euclidean.dist))
system.time(Dy <- symm(X, euclidean.dist, symm = TRUE))
system.time(Dv <- w_vapp(X, euclidean.dist))
