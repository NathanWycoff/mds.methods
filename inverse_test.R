source('forward_methods.R')
source("inverse_methods.R")

##Generate some weights and try to recover them.
#Some params
n <- 20#How many points?
p <- 6#How many high D dims?
k <- 2#How many low D dims (never not 2, but I like to be overly general)
n.inits <- 1
dist.func <- inner.prod.dist
seed <- 123

#Generate the weights and induced lowD points
high_d <- matrix(rnorm(n*p), ncol = p)
user_weights <- rgamma(p, 1, 1)
user_weights <- user_weights / sum(user_weights)
low_d_result <- forward_mds(high_d, k, user_weights, dist.func, n.inits, seed + 1)
user_low_d <- low_d_result$par

#Try to recover those weights based on the low D points
forward.n.inits <- n.inits#What we used to call n.inits changes
n.inits <- 2#How many times do we try new weights?
inferred_weights <- inverse_step(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed) 
