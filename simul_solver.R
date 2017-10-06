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

approx_smacof_inverse_step <- function(user_low_d, prev_low_d, high_d, prev_weights, dist.func, forward.n.inits, seed) {
    set.seed(seed)

    #Calculate the low d distance matrix
    user_low_d_dist <- good.dist(user_low_d, dist.func)

    #Get some random
    p <- dim(high_d)[2]

    #Run the solver a couple times
    result <- optim(prev_weights, approx_smacof_inverse_cost, method = "L-BFGS-B",
                    user_low_d_dist = user_low_d_dist, high_d = high_d, 
                    prev_low_d = prev_low_d, n.inits = forward.n.inits, 
                    dist.func = dist.func, lower = 0.01)

    return(result)
}
