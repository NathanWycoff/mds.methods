#!/usr/bin/Rscript
#  do_we_match_smacof.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.29.2017

##This script tests whether we match the results of the smacof package from CRAN.


#######Import the source code
source("forward_methods.R")

###Do mds from the SMACOF package with multiple random initializations
mult.init.mds <- function(n.inits, ...) {
    sols <- lapply(1:n.inits, function(i) mds(...))
    optim.sol <- sols[[which.min(sapply(sols, function(i) i$stress))]]
    return(optim.sol)
}

########Create some data
seed <- 123
set.seed(seed)
p <- 30
n <- 120
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL

######## Compare on arbitrary points
require(smacof)
seed <- 123
weights = rep(1,p)
system.time(low_d_us <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist))

delta <- good.dist(high_d, euclidean.dist) 

system.time(low_d_them <- mult.init.mds(n.inits = 100, delta, init = 'random'))

#Get and compare distance matrices
us_dist <- as.matrix(dist(low_d_us$par))
them_dist <- as.matrix(dist(low_d_them$conf))
us_dist <- us_dist / sum(us_dist)
them_dist <- them_dist / sum(them_dist)
print(norm(us_dist - them_dist, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_d_us$par, main=low_d_us$value, lty=0, lwd=0)
text(low_d_us$par[,1], low_d_us$par[,2], 1:n)
plot(low_d_them$conf, main=low_d_them$stress, lty=0, lwd=0)
text(low_d_them$conf[,1], low_d_them$conf[,2], 1:n)


######## Next, see if we can match PCA using inner product distances.
## This will test if the smacof package's function takes diagonal elements of the
# distance matrix into account.
require(smacof)
seed <- 123
weights = rep(1,p)

to_plot <- data.frame(prcomp(high_d)$x[,1:2])

delta <- good.dist(high_d, inner.prod.dist) 

system.time(low_d_them <- mult.init.mds(n.inits = 100, delta, init = 'random'))

#Get and compare distance matrices
us_dist <- as.matrix(dist(low_d_pca$par))
them_dist <- as.matrix(dist(low_d_them$conf))
us_dist <- us_dist / sum(us_dist)
them_dist <- them_dist / sum(them_dist)
print(norm(us_dist - them_dist, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_d_pca$par, main=low_d_us$value, lty=0, lwd=0)
text(low_d_pca$par[,1], low_d_us$par[,2], 1:n)
plot(low_d_them$conf, main=low_d_them$stress, lty=0, lwd=0)
text(low_d_them$conf[,1], low_d_them$conf[,2], 1:n)
