#!/usr/bin/Rscript
#  forward_methods.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.29.2017

#TODO: Make good.distance accept weights; remove all instances of weighting outside of the function.

#################################################################################
##### Functions for forwards
#An inner product function
inner.prod.dist <- function(x, y) t(x) %*% y
euclidean.dist <- function(x, y) sqrt(sum((x-y)^2))
manhattan.dist <- function(x, y) sum(abs(x-y))

#A distance function that:
#1) Returns a matrix instead of a dist object (I need the diagonals for inner prod distance)
#2) Allows for custom distance function specification.
good.dist <- function(X, dist.func, weights = NULL) {
    if (!is.null(weights)) {
        weights <- sqrt(weights)
        X <- X %*% diag(weights)
    }
    n <- dim(X)[1]
    dists <- matrix(0, ncol = n, nrow = n)
    for (i in 1:n) {
        for (j in 1:n) {
            dists[i,j] <- dist.func(X[i,], X[j,])
        }
    }
    return(dists)
}

# Can standardize the stress by dividing by the sum of high dim distances squared.
# Argument is std that takes on values TRUE or FALSE.

forward_cost <- function(low_d, high_d_dist, dist.func, std) {
  low_d <- matrix(low_d, ncol = 2)
  low_d_dist <- good.dist(low_d, dist.func)
  
  #Normalize the distance matrices
  high_d_dist <- high_d_dist / sum(high_d_dist)
  low_d_dist <- low_d_dist / sum(low_d_dist)
  
  diff_mat <- low_d_dist - high_d_dist
  
  #Calculate the stress of all elements
  stress <- sum(diff_mat^2)
  
  ## Standardize the stress
  stress_s <- stress/sum(high_d_dist^2)
  
  if(std==TRUE){return(stress_s)}
  else{return(stress)}
}

forward_mds  <- function(high_d, weights, dist.func, 
                 thresh = 1e-5, max.iters = 1000, n.inits = 10, 
                 seed = NULL, std=TRUE) {
  if(is.null(seed)) {
    seed <- sample(1:100000, 1)
  }
  
  #Create the distance matrix
  n <- dim(high_d)[1]
  high_d_dist <- good.dist(high_d, dist.func, weights)
  
  #Randomly Init Points
  set.seed(seed)
  inits <- lapply(1:n.inits, function(x) rnorm(n*2))
  
  #Get optimal points
  z_optimals <- lapply(inits, function(init) optim(init, forward_cost, method = "BFGS", high_d_dist = high_d_dist, dist.func = dist.func, std=std, control = list('abstol' = 1e-15)))
  
  costs <- unlist(lapply(z_optimals, function(z) z$value))
  optimal <- which.min(costs)
  z_optimal <- z_optimals[[optimal]]
  
  #Shape it into a matrix
  z_optimal$par <- matrix(z_optimal$par, ncol = 2)
  z_optimal$par <- scale(z_optimal$par)
  attr(z_optimal$par, 'scaled:center') <- NULL
  attr(z_optimal$par, 'scaled:scale') <- NULL
  
  return(z_optimal)
}


##################################################SMACOF funcs
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

smacof_forward_mds <- function(high_d, weights, dist.func = euclidean.dist,
                  thresh = 1e-5, max.iters = 1000, n.inits = 10, seed = NULL, std=TRUE) {
  if (is.null(seed)) {
    seed <- sample(1:10000, 1)
  }
  set.seed(seed)
  
  #Get our high D distance
  if (sum(weights) < 1e-5) {
    stop("At least some weights should be nonnegative")
  }
  true_dist <- good.dist(high_d, dist.func, weights)
  true_dist <- true_dist / sum(true_dist)
  
  #Run a bunch of smacof algos
  results <- lapply(1:n.inits, function(i)
    single_smacof(true_dist = true_dist, dist.func = dist.func,
                  thresh = thresh, max.iters = max.iters, std=std))
  
  #Get lowest cost result
  costs <- sapply(results, function(i) i$value)
  
  return(results[[which.min(costs)]])
}


single_smacof <- function(true_dist, dist.func = euclidean.dist,
                          thresh = 1e-5, max.iters = 1000, init_low_d = NULL, std=TRUE) {
  #Make a random initial lowD matrix if non is provided
  if (is.null(init_low_d)) {
    n <- dim(true_dist)[1]
    low_d <- matrix(rnorm(n*2), ncol = 2)
  } else{
    low_d <- init_low_d
  }
  
  ##Do a SMACOF iteration
  diff <- Inf
  last_err <- Inf
  iter <- 0
  while (diff > thresh && iter < max.iters) {
    iter <- iter + 1
    
    #Create the B matrix
    B <- make.B(true_dist, low_d)
    
    #Get our next Z
    low_d <- B %*% low_d
    
    #Print current error
    err <- forward_cost(low_d, true_dist, dist.func, std=std)
    
    #Get difference in error
    diff <- abs(last_err - err)
    last_err <- err
  }
  
  if (iter == max.iters) {
    warning("SMAMOF algo did not converge")
  }
  
  #Scale the coords
  low_d <- scale(low_d)
  attr(low_d, 'scaled:center') <- NULL
  attr(low_d, 'scaled:scale') <- NULL
  err <- forward_cost(low_d, true_dist, dist.func, std=std)
  
  return(list(par = low_d, value = err))
}
