#!/usr/bin/Rscript
#  inverse_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.03.2017
source('forward_methods.R')
source("inverse_methods.R")

##Generate some weights and try to recover them.
#Some params
#Create some data
seed <- 123
set.seed(seed)
p <- 5
n <- 5
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL
k <- 2#How many low D dims (never not 2, but I like to be overly general)

#Specify some params
n.inits <- 10
dist.func <- euclidean.dist

#Generate the weights and induced lowD points
user_weights <- rgamma(p, 1, 1)*10
user_weights <- user_weights / sum(user_weights)
low_d_result <- forward_mds(high_d, k, user_weights, dist.func, n.inits, seed + 1)
user_low_d <- low_d_result$par

#Try to recover those weights based on the low D points
set.seed(seed)
forward.n.inits <- n.inits#What we used to call n.inits changes
system.time(inferred_weights <- inverse_step(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed))

print(inferred_weights$par - user_weights)

#Now try to do it using the SMACOF method on the forward, still numerically tho for the backward
#Try to recover those weights based on the low D points
set.seed(seed)
forward.n.inits <- n.inits#What we used to call n.inits changes
system.time(inferred_weights <- smacof_inverse_mds(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed + 1))

print(inferred_weights$par - user_weights)
