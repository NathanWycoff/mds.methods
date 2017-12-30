#!/usr/bin/Rscript
#  simul_solver_again.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.09.2017
source("forward_methods.R")
source("inverse_methods.R")

w_step <- function(low_d, last_low_d, high_d, dist.func, prev_weights) {
    #Make low D dist, and normalize it
    low_d_dist <- good.dist(low_d, dist.func)
    low_d_dist <- low_d_dist / sum(low_d_dist)

    #Scale the high D data
    high_d <- scale(high_d)
    attr(high_d, 'scaled:center') <- NULL
    attr(high_d, 'scaled:scale') <- NULL

    #Create our cost function
    high_d_dist <- function(weights) {
        high_d_dist <- good.dist(high_d, dist.func, weights)
        high_d_dist <- high_d_dist / sum(high_d_dist)
        return(high_d_dist)
    }

    #Get the approximate solution based on weights (one iter of SMACOF)
    approx_sol <- function(weights) {
        Z <- make.B(high_d_dist(weights), last_low_d) %*% last_low_d
        Z <- scale(Z)
        attr(Z, 'scaled:center') <- NULL
        attr(Z, 'scaled:scale') <- NULL
        return(Z)
    }

    #Put everything together in our cost function
    cost <- function(weights) {
        induced_dist <- good.dist(approx_sol(weights), dist.func)
        induced_dist <- induced_dist / sum(induced_dist)
        return(sum((low_d_dist - induced_dist)^2))
    }
    
    #Do the optimization numerically, for now.
    optim_weights <- optim(prev_weights, cost, method = 'L-BFGS-B', 
                           lower = 0.00001)

    return(optim_weights)
}


w_single_step <- function(low_d, last_low_d, high_d, dist.func, prev_weights) {
    #Make low D dist, and normalize it
    low_d_dist <- good.dist(low_d, dist.func)
    low_d_dist <- low_d_dist / sum(low_d_dist)

    #Scale the high D data
    high_d <- scale(high_d)
    attr(high_d, 'scaled:center') <- NULL
    attr(high_d, 'scaled:scale') <- NULL

    #Create our cost function
    high_d_dist <- function(weights) {
        high_d_dist <- good.dist(high_d, dist.func, weights)
        high_d_dist <- high_d_dist / sum(high_d_dist)
        return(high_d_dist)
    }

    #Get the approximate solution based on weights (one iter of SMACOF)
    approx_sol <- function(weights) {
        Z <- make.B(high_d_dist(weights), last_low_d) %*% last_low_d
        Z <- scale(Z)
        attr(Z, 'scaled:center') <- NULL
        attr(Z, 'scaled:scale') <- NULL
        return(Z)
    }

    #Put everything together in our cost function
    cost <- function(weights) {
        induced_dist <- good.dist(approx_sol(weights), dist.func)
        induced_dist <- induced_dist / sum(induced_dist)
        return(sum((low_d_dist - induced_dist)^2))
    }
    
    #Do the optimization numerically, for now.
    optim_weights <- optim(prev_weights, cost, method = 'L-BFGS-B', 
                           lower = 0.00001, control = list('maxit'=1))

    return(optim_weights)
}

#Do a SMACOF step to update Z
Z_step <- function(prev_low_d, current_weights, high_d) {
    #Get high D distance
    high_d_dist <- good.dist(high_d, dist.func, current_weights)
    
    new_low_d <- make.B(high_d_dist, prev_low_d) %*% prev_low_d

    return(new_low_d)
}

#Iterate between the two
approx_inverse_smacof <- function(low_d, last_low_d, high_d, init_weights,
                                dist.func, thresh = 1e-5, max.iters = 1e4) {
    #Scale the high D data
    high_d <- scale(high_d)
    attr(high_d, 'scaled:center') <- NULL
    attr(high_d, 'scaled:scale') <- NULL

    #Initialize things
    weights <- init_weights
    new_weights <- weights
    Z <- last_low_d

    diff <- Inf
    iter <- 0
    while (diff > thresh && iter < max.iters) {
        iter <- iter + 1
        #Store the old weights
        old_weights <- new_weights

        #Update weights
        lambda <- 0.9
        new_weights <- w_single_step(low_d, Z, high_d, dist.func, old_weights)$par
        new_weights <- new_weights / sum(new_weights)
        weights <- lambda*weights + (1-lambda) * new_weights
        weights <- weights / sum(weights)

        #Update low D points
        new_Z <- Z_step(Z, weights, high_d)
        new_Z <- scale(new_Z)
        Z <- lambda * Z + (1-lambda) * new_Z
        #high_d_dist <- good.dist(high_d, dist.func, weights)
        #Z <- single_smacof(high_d_dist, dist.func=dist.func, init_low_d=Z)$par

        #See if we're converged
        diff <- norm(new_weights - old_weights, '2')
    }

    if (iter == max.iters) {
        warning("approx inverse smacof did not converge")
    }

    return(list(Z = Z, weights = weights))
}
