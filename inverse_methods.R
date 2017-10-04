#!/usr/bin/Rscript
#  inverse_methods.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.03.2017
source("forward_methods.R")

################################################################################
## Functions for inverse
inverse_cost <- function(weights, user_low_d_dist, high_d, k, n.inits, dist.func) {
    #Get the low d projection induced by the weights, and its distance matrix
    weights_low_d <- forward_mds(high_d, k, weights, dist.func, n.inits, sample(1:1000,1))$par
    weights_low_d_dist <- good.dist(weights_low_d, dist.func)

    #Compare the two distance matrices
    diff <- weights_low_d_dist - user_low_d_dist
    stress <- sum(diff^2)

    return(stress)
}

###
#@param n.inits the number of intial weight configs
#@param forward.n.inits the number of initial point configs for each forward step
inverse_step <- function(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed) {
    set.seed(seed)

    #Calculate the low d distance matrix
    user_low_d_dist <- good.dist(user_low_d, dist.func)

    #Get some random
    p <- dim(high_d)[2]
    init.weights <- lapply(1:n.inits, function(i) rgamma(p, 1, 1))

    #Run the solver a couple times
    results <- lapply(init.weights, function(init) optim(init, inverse_cost, method = "L-BFGS-B", user_low_d_dist = user_low_d_dist, high_d = high_d, k = k, n.inits = forward.n.inits, dist.func = dist.func, lower = 0))
    result <- results[[which.min(sapply(results, function(i) i$value))]]

    return(result)
}


########################################################## SMACOF methods
smacof_inverse_cost <- function(weights, user_low_d_dist, high_d, k, n.inits, dist.func) {
    #Get the low d projection induced by the weights, and its distance matrix
    weights_low_d <- forward_mds(high_d, k, weights, dist.func, n.inits, sample(1:1000,1))$par
    weights_low_d_dist <- good.dist(weights_low_d, dist.func)

    #Compare the two distance matrices
    diff <- weights_low_d_dist - user_low_d_dist
    stress <- sum(diff^2)

    return(stress)
}

smacof_inverse_mds <- function(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed) {
    set.seed(seed)

    #Calculate the low d distance matrix
    user_low_d_dist <- good.dist(user_low_d, dist.func)

    #Get some random
    p <- dim(high_d)[2]
    init.weights <- lapply(1:n.inits, function(i) rgamma(p, 1, 1))

    #Run the solver a couple times
    results <- lapply(init.weights, function(init) optim(init, smacof_inverse_cost, method = "L-BFGS-B", user_low_d_dist = user_low_d_dist, high_d = high_d, k = k, n.inits = forward.n.inits, dist.func = dist.func, lower = 0))
    result <- results[[which.min(sapply(results, function(i) i$value))]]

    return(result)
}

