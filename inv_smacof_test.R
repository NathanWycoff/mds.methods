source("smacof.R")

#Stability under perturbations of the weight vector
orig <- smacof_forward_mds(X, weights = weights, n.inits = 100, dist.func = euclidean.dist)
weights <- weights + rnorm(p,0,0.001)
perturbed <- smacof_forward_mds(X, weights = weights, n.inits = 100, dist.func = euclidean.dist)

par(mfrow=c(2,1))
plot(orig$low_d)
plot(perturbed$low_d)

#Generate some low_d points
set.seed(123)
weights <- rgamma(p, 1, 1)
user_low_d <- smacof_forward_mds(X, weights = weights, n.inits = 100, dist.func = euclidean.dist)$low_d


smacof_inverse_cost <- function(weights, user_low_d_dist, high_d, k = 2, n.inits, dist.func) {
    #Get the low d projection induced by the weights, and its distance matrix
    weights_low_d <- smacof_forward_mds(high_d, weights, dist.func, n.inits = n.inits)$low_d
    weights_low_d_dist <- good.dist(weights_low_d, dist.func)

    #Compare the two distance matrices
    diff <- weights_low_d_dist - user_low_d_dist
    stress <- sum(diff^2)

    return(stress)
}

###
#@param n.inits the number of intial weight configs
#@param forward.n.inits the number of initial point configs for each forward step
smacof_inverse_step <- function(user_low_d, high_d, dist.func, n.inits, forward.n.inits, seed) {
    set.seed(seed)

    #Calculate the low d distance matrix
    user_low_d_dist <- good.dist(user_low_d, dist.func)

    #Get some random
    p <- dim(high_d)[2]
    init.weights <- lapply(1:n.inits, function(i) rgamma(p, 1, 1))

    #Run the solver a couple times
    init <- init.weights[[1]] 
    result <- optim(init, smacof_inverse_cost, method = "L-BFGS-B", user_low_d_dist = user_low_d_dist, high_d = high_d, k = k, n.inits = forward.n.inits, dist.func = dist.func, lower = 0)

    return(result)
}
