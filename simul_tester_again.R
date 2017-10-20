#!/usr/bin/Rscript
#  inverse_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.03.2017
require(MCMCpack)
source('forward_methods.R')
source("inverse_methods.R")
source('simul_solver_again.R')

##Generate some weights and try to recover them.
#Some params
thresh <- 1e-5
max.iters <- 1000

#Create some data
seed <- 1234
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


############################# Test that starting with the right solution stays there
#Generate the 'initial' data
set.seed(seed)
init_weights <- rgamma(p, 1, 1)
init_weights <- init_weights / sum(init_weights)
low_d_result <- smacof_forward_mds(high_d, init_weights, dist.func, n.inits, seed + 1)
init_low_d <- low_d_result$par

####Do the actual recovering
w_step(init_low_d, init_low_d, high_d, dist.func, init_weights)

############################# Try to get a slight perturbation of the weights
#Generate the 'initial' data
set.seed(seed)
init_weights <- rgamma(p, 1, 1)
init_weights <- init_weights / sum(init_weights)
low_d_result <- smacof_forward_mds(high_d, init_weights, dist.func, n.inits, seed + 1)
init_low_d <- low_d_result$par

#Different weights
#perturb_weights <- rep(1,p)
#perturb_weights[2] <- 10
#perturb_weights <- rgamma(p,1,1)
set.seed(seed+1)
perturb_weights <- init_weights * rgamma(p, 10, 10)
perturb_weights <- perturb_weights / sum(perturb_weights)
low_d_result <- smacof_forward_mds(high_d, perturb_weights, dist.func, n.inits, seed + 1)
perturb_low_d <- low_d_result$par

#Plot the two results
par(mfrow=c(2,1))
plot(init_low_d, main='init', lty=0, lwd=0)
text(init_low_d[,1], init_low_d[,2], 1:n)
plot(perturb_low_d, main='perturb', lty=0, lwd=0)
text(perturb_low_d[,1], perturb_low_d[,2], 1:n)

####Do the actual recovering
result <- w_step(perturb_low_d, init_low_d, high_d, dist.func, init_weights)
infer_weights <- result$par / sum(result$par)

#Does it do any better than the initial solution?
print(infer_weights)
print(perturb_weights)
print(init_weights)

print(sum((infer_weights - perturb_weights)^2))
print(sum((init_weights - perturb_weights)^2))

############################# Try to get entirely new weights using the approx
#Generate the 'initial' data
set.seed(seed)

mine <- c()
old <- c()
rand <- c()

iters <- 100

for (i in 1:iters) {
    alpha <- 0.1#Dirichlet concentration parameter
    weights_1 <- as.numeric(rdirichlet(1, rep(alpha, p)))
    weights_1 <- weights_1 / sum(weights_1)
    low_d_result <- smacof_forward_mds(high_d, weights_1, dist.func, n.inits, seed + 1)
    low_d_1 <- low_d_result$par

    #Different weights
    #weights_2 <- rep(1,p)
    #weights_2[2] <- 10
    weights_2 <- as.numeric(rdirichlet(1, rep(alpha, p)))
    weights_2 <- weights_2 / sum(weights_2)
    low_d_result <- smacof_forward_mds(high_d, weights_2, dist.func, n.inits, seed + 1)
    low_d_2 <- low_d_result$par

    #Plot the two results
    par(mfrow=c(2,1))
    plot(low_d_1, main='first', lty=0, lwd=0)
    text(low_d_1[,1], low_d_1[,2], 1:n)
    plot(low_d_2, main='second', lty=0, lwd=0)
    text(low_d_2[,1], low_d_2[,2], 1:n)

    ####Do the actual recovering
    result <- approx_inverse_smacof(low_d_2, low_d_1, high_d, weights_1,
                                    dist.func = dist.func)
    infer_weights <- result$weights / sum(result$weights)

    old_result <- old_inverse(low_d_1, high_d, dist.func)
    old_weights <- old_result$par / sum(old_result$par)


    #Does it do any better than the initial solution?
    mine <- c(mine, mean(abs(infer_weights - weights_2)))
    old <- c(old, mean(abs(old_weights - weights_2)))
    rand <- c(rand, mean(abs(weights_1 - weights_2)))

    #print(sum(abs(infer_weights - weights_2)))
    #print(sum(abs(old_weights - weights_2)))
    #print(sum(abs(weights_1 - weights_2)))
    ##See if the inferred Z's look anything like the User's
    #par(mfrow=c(2,2))
    #plot(result$Z, main='Estimated', lty=0, lwd=0)
    #text(result$Z[,1], result$Z[,2], 1:n)
    #plot(low_d_2, main='User', lty=0, lwd=0)
    #text(low_d_2[,1], low_d_2[,2], 1:n)
    #plot(low_d_1, main='Original', lty=0, lwd=0)
    #text(low_d_1[,1], low_d_1[,2], 1:n)
}

print(mean(mine))
print(mean(old))
print(mean(rand))
