#!/usr/bin/Rscript
#  simul_solver.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.04.2017

#Creates the transform from the previous point
make.B <- function(true_dist, low_d) {
    #Create the low D distance matrix
    low_dist <- as.matrix(dist(low_d))

    n <- dim(true_dist)[1]
    B <- matrix(0, ncol = n, nrow = n)
    
    #iterate over the B matrix to create it
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                B[i,j] <- -true_dist[i,j] / low_dist[i,j]
            }
        }
    }

    #Create the diagonal elements
    diag(B) <- -colSums(B)

    return(B)
}

##Get the cost for the approximate SMACOF-based inverse MDS algorithm
approx_smacof_inverse_cost <- function(weights, user_low_d_dist, 
                                 high_d, prev_low_d, n.inits, dist.func) {
    #Get the low d projection induced by the weights, and its distance matrix
    #This next line is the one that needs to change
    #weights_low_d <- forward_mds(high_d, k, weights, dist.func, n.inits, sample(1:1000,1))$par
    #Get the high d distance
    weights <- sqrt(weights)
    high_d <- high_d %*% diag(weights)
    true_dist <- good.dist(high_d, dist.func)

    #Get the approximate solution
    weights_low_d <- 1/n * make.B(true_dist, prev_low_d) %*% prev_low_d
    weights_low_d_dist <- good.dist(weights_low_d, dist.func)

    #Compare the two distance matrices
    diff <- weights_low_d_dist - user_low_d_dist
    stress <- sum(diff^2)

    return(stress)
}

#Opmitize the weights for 1 step of the approximate SMACOF algo
approx_smacof_inverse_step <- function(user_low_d, prev_low_d, high_d, 
                                       prev_weights, dist.func, forward.n.inits, 
                                       seed) {
    set.seed(seed)

    #Calculate the low d distance matrix
    user_low_d_dist <- good.dist(user_low_d, dist.func)

    p <- dim(high_d)[2]

    #Run the solver a couple times
    result <- optim(prev_weights, approx_smacof_inverse_cost, method = "L-BFGS-B",
                    user_low_d_dist = user_low_d_dist, high_d = high_d, 
                    prev_low_d = prev_low_d, n.inits = forward.n.inits, 
                    dist.func = dist.func, lower = 0.0001)

    return(result)
}

approx_smacof_simul <- function(user_low_d, prev_low_d, high_d, 
                               prev_weights, dist.func, forward.n.inits, thresh,
                               max.iters, seed) {

    #Initialize on the last values
    weights <- prev_weights
    Z <- prev_low_d

    #Iterate between optimizing weights and low_d coords
    diff <- Inf
    last_err <- Inf
    iter <- 0
    while (diff > thresh && iter < max.iters) {
        in_seed <- sample(1:10000, 1)

        #Update weights
        last_weights <- weights
        weights <- approx_smacof_inverse_step(user_low_d, Z, high_d, weights,
                                              dist.func, forward.n.inits,
                                              in_seed)$par
        weights <- weights / sum(weights)
        
        #Update low D points
        w_high_d <- high_d %*% diag(sqrt(weights))
        true_dist <- good.dist(w_high_d, dist.func)
        Z <- single_smacof(true_dist, dist.func = dist.func, thresh = thresh, 
                           max.iters = max.iters)$par

        #Update iterative params
        diff <- norm(weights - last_weights, '2')
        iter <- iter + 1
    }

    return(list(weights = weights, Z = Z))
}
