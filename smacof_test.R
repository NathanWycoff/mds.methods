#Import the source code
source("forward_methods.R")
source("smacof.R")

#Create some data
set.seed(123)
p <- 30
n <- 30
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL

#Compare two runs of smacof
weights <- rep(1, p)
low_1 <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist)
weights <- rep(1, p)
low_2 <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist)

#Get and compare distance matrices
smacof_dist <- as.matrix(dist(low_1$low_d))
optim_dist <- as.matrix(dist(low_2$low_d))
print(norm(smacof_dist - optim_dist, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_1$low_d, main=low_1$cost, lty=0, lwd=0)
text(low_1$low_d[,1], low_1$low_d[,2], 1:n)
plot(low_2$low_d, main=low_2$cost, lty=0, lwd=0)
text(low_2$low_d[,1], low_2$low_d[,2], 1:n)

#Check that smacof and optim give the same results with uniform weights
seed <- 123
weights = rep(1, p)
system.time(low_d_smacof <- smacof_forward_mds(high_d, weights = weights, 
                                n.inits = 100, dist.func = euclidean.dist))
system.time(low_d_optim <- forward_mds(high_d, k = 2, weights = weights, 
                                dist.func = euclidean.dist, n.inits = 100, seed))

#Get and compare distance matrices
smacof_dist <- as.matrix(dist(low_d_smacof$low_d))
optim_dist <- as.matrix(dist(low_d_optim$par))
print(norm(smacof_dist - optim_dist, '2'))

#Plot the results
par(mfrow=c(2,1))
plot(low_d_smacof$low_d, main=low_d_smacof$cost, lty=0, lwd=0)
text(low_d_smacof$low_d[,1], low_d_smacof$low_d[,2], 1:n)
plot(low_d_optim$par, main=low_d_optim$value, lty=0, lwd=0)
text(low_d_optim$par[,1], low_d_optim$par[,2], 1:n)
