#!/usr/bin/Rscript
#  smacof_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.03.2017

#Import the source code
source("../R/forward_methods.R")

#Create some data
seed <- 123
set.seed(seed)
p <- 30
n <- 6
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL

######################################## Dumb MDS
n.inits <- 10
dist.func <- euclidean.dist

#Try the algo with two different seeds and see if it converges to the same thing
weights <- rep(1,p)
low_d_result <- forward_mds(high_d, weights, dist.func, n.inits = n.inits, 
                            seed = seed)
low_d <- low_d_result$par

par(mfrow=c(2,1))
plot(low_d, lwd = 0, main = paste('MDS baby', 'stress:',low_d_result$value), lty=0)
text(low_d[,1], low_d[,2], 1:n)

old_low_d <- low_d

low_d_result <- forward_mds(high_d, weights, dist.func, n.inits = n.inits,
                           seed =  seed + 1)
low_d <- low_d_result$par

plot(low_d, lwd = 0, main = paste('MDS baby', 'stress:',low_d_result$value), lty=0)
text(low_d[,1], low_d[,2], 1:n)

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d))
low_2 <- as.matrix(dist(old_low_d))
print(norm(low_1 - low_2))

###################################### SMACOF MDS
#Compare two runs of smacof
weights <- rep(1, p)
low_1 <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist)
weights <- rep(1, p)
low_2 <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist)

#Get and compare distance matrices
dist_1 <- as.matrix(dist(low_1$par))
dist_2 <- as.matrix(dist(low_2$par))
print(norm(dist_1 - dist_2, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_1$par, main=low_1$value, lty=0, lwd=0)
text(low_1$par[,1], low_1$par[,2], 1:n)
plot(low_2$par, main=low_2$value, lty=0, lwd=0)
text(low_2$par[,1], low_2$par[,2], 1:n)

###################################### Comparison
#Check that smacof and optim give the same results with random weights
seed <- 123
weights = rgamma(p,1,1)
system.time(low_d_smacof <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist))
system.time(low_d_optim <- forward_mds(high_d, weights = weights, 
                                dist.func = euclidean.dist, n.inits = 100, 
                                seed = seed))

#Get and compare distance matrices
smacof_dist <- as.matrix(dist(low_d_smacof$par))
optim_dist <- as.matrix(dist(low_d_optim$par))
smacof_dist <- smacof_dist / sum(smacof_dist)
optim_dist <- optim_dist / sum(optim_dist)
print(norm(smacof_dist - optim_dist, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_d_smacof$par, main=low_d_smacof$value, lty=0, lwd=0)
text(low_d_smacof$par[,1], low_d_smacof$par[,2], 1:n)
plot(low_d_optim$par, main=low_d_optim$value, lty=0, lwd=0)
text(low_d_optim$par[,1], low_d_optim$par[,2], 1:n)
