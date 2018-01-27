#!/usr/bin/Rscript
#  forward_methods.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.29.2017

#TODO: We should separate low D high D distances.
#TODO: We should use a stopping threshold that is relative as opposed to absolute.

#################################################################################
##### Functions for forwards
#An inner product function
inner.prod.dist <- function(x, y) t(x) %*% y
euclidean.dist <- function(x, y) sqrt(sum((x-y)^2))
manhattan.dist <- function(x, y) sum(abs(x-y))

#' A distance function that, unlike the built-in dist:
#' 1) Returns a matrix instead of a dist object (I need the diagonals for inner prod distance)
#' 2) Allows for custom distance function specification.
#' @param X An n by p real valued matrix, rows representing observations, columns representing features. A distance matrix will be computed along its rows.
#' @param dist.func A function which accepts two arguments: two rows X, and returns a distance (or generalized distance) between them.
#' @param weights Either a length p nonnegative vector or NULL. If not null, the ith column of X will be scaled by the root of the ith element of weights prior to distance calculation. If not, no scaling will occur.
#' @param symm Is the distance function symmetric? If so, we can skip a bunch of function evaluations. Default is FALSE.
#' @return An n by n matrix of distances, element i,j of which being the result of dist.func applied to rows i and j.
#' @export
#' @examples 
good.dist <- function(X, dist.func, weights = NULL, symm = FALSE) {
    if (!is.null(weights)) {
        weights <- sqrt(weights)
        X <- X %*% diag(weights)
    }
    n <- dim(X)[1]
    dists <- matrix(0, ncol = n, nrow = n)

    #Calculate the actual distances
    if (symm) {
        for (i in 1:n) {
            for (j in 1:i) {
                dists[i,j] <- dists[j,i] <- dist.func(X[i,], X[j,])
            }
        } 
    } else {
        for (i in 1:n) {
            for (j in 1:n) {
                dists[i,j] <- dist.func(X[i,], X[j,])
            }
        }
    }
    return(dists)
}

#' Returns MDS stress for the forward algorithm as defined in equation 8.15 of Borg and Groenen, with all weights (w_{i,j}) equal to 1 (do not confuse these with the weights on dimensions defined in this package).
#' @param low_d The low_d solution, an n by 2 matrix, the cost of which is to be evaluated.
#' @param high_d_dist The distance matrix, n by n, of the high dimensional data of which we seek a low D approximation.
#' @param std Boolean, should stress be standardized? This is accomplished by dividing stress by the sum of squared high D distance
#' @return Nonnegative scalar stress of the configuration
#' @export
#' @examples 
forward_cost <- function(low_d, high_d_dist, std = TRUE) {
    low_d <- matrix(low_d, ncol = 2)
    low_d_dist <- good.dist(low_d, euclidean.dist, symm = TRUE)

    #Normalize the distance matrices
    high_d_dist <- high_d_dist / sum(high_d_dist)
    low_d_dist <- low_d_dist / sum(low_d_dist)


    diff_mat <- low_d_dist - high_d_dist

    #Calculate the stress of all elements
    stress <- sum(diff_mat^2)

    #Standardize if necessary
    if (std)  {
        stress <- stress/sum(high_d_dist^2)
    }

    return(stress)
}


#' Do forward MDS by numerically minimizing the stress defined in forward_cost.
#' @param high_d The high dimensional data of which a low dimensional representation is desired, an n by p matrix where rows represent observations
#' @param low_d The low_d solution, an n by 2 matrix, the cost of which is to be evaluated.
#' @param dist.func The distance function to be used for high D distance computation; euclidean is always used for low D distance.
#' @param thresh The threshold below which stress must fall before computation ends. 
#' @param max.iters The maximum number of iterations passed to optim for stress minimization.
#' @param n.inits The number of times the optim is run from a random configuration. This can be important, as the cost surface is highly nonconvex.
#' @param seed Random seed used for initialization
#' @param std Boolean, should stress be standardized? This is accomplished by dividing stress by the sum of squared high D distance
#' @param symm Boolean, is the distance function symmetric? We can save on computation if so.
#' @return List with $par, the optimal configuration as an n by 2 matrix, and $value as the stress of this configuration
#' @export
#' @examples 
forward_mds <- function(high_d, weights, dist.func, 
                   thresh = 1e-5, max.iters = 1000, n.inits = 10, 
                   seed = NULL, std = TRUE, symm = FALSE) {
    if(is.null(seed)) {
        seed <- sample(1:100000, 1)
    }

    #Create the distance matrix
    n <- dim(high_d)[1]
    high_d_dist <- good.dist(high_d, dist.func, weights, symm)
    
    #Randomly Init Points
    set.seed(seed)
    inits <- lapply(1:n.inits, function(x) rnorm(n*2))

    #Get candidate solutions
    z_optimals <- lapply(inits, function(init) optim(init, forward_cost, method = "BFGS", 
        std = std, high_d_dist = high_d_dist, 
        control = list('abstol' = thresh, 'maxit' = max.iters)))

    #Find the optimal sol among the candidates
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

#' Do forward MDS by using the SMACOF algorithm defined on page 191 of Borg and Groenen.
#' @param high_d The high dimensional data of which a low dimensional representation is desired, an n by p matrix where rows represent observations
#' @param low_d The low_d solution, an n by 2 matrix, the cost of which is to be evaluated.
#' @param dist.func The distance function to be used for both low and high D distance computation.
#' @param thresh The threshold below which stress must fall before computation ends. 
#' @param max.iters The maximum number of SMACOF iterations (i.e., the maximum number of times we solve the quadratic system described in chapter 8 of Borg and Groenen) per initialization.
#' @param n.inits The number of times the SMACOF algorithm is run from a random configuration. This can be important, as the cost surface is highly nonconvex.
#' @param seed Random seed used for initialization
#' @param symm Boolean, is the distance function symmetric? We can save on computation if so.
#' @return List with $par, the optimal configuration as an n by 2 matrix, and $value as the stress of this configuration
#' @export
smacof_forward_mds <- function(high_d, weights, dist.func = euclidean.dist,
                   thresh = 1e-5, max.iters = 1000, n.inits = 10, seed = NULL,
                   std = TRUE, symm = FALSE) {
    if (is.null(seed)) {
        seed <- sample(1:10000, 1)
    }
    set.seed(seed)

    #Get our high D distance
    if (sum(weights) < 1e-5) {
        stop("At least some weights should be actually positive; 
             their sum seems to be very small")
    }
    true_dist <- good.dist(high_d, dist.func, weights, symm)
    true_dist <- true_dist / sum(true_dist)

    #Run a bunch of smacof algos
    results <- lapply(1:n.inits, function(i)
                    single_smacof(true_dist = true_dist, dist.func = dist.func,
                        thresh = thresh, max.iters = max.iters))

    #Get lowest cost result
    costs <- sapply(results, function(i) i$value)

    return(results[[which.min(costs)]])
}

single_smacof <- function(true_dist, dist.func = euclidean.dist,
                   thresh = 1e-5, max.iters = 1000, init_low_d = NULL,
                   std = TRUE) {
    #Make a random initial lowD matrix if none is provided
    n <- nrow(true_dist)
    if (is.null(init_low_d)) {
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
        err <- forward_cost(low_d, true_dist, std = std)

        #Get difference in error
        diff <- abs(last_err - err)
        last_err <- err
    }

    if (iter == max.iters) {
        warning("SMACOF algo did not converge")
    }

    #Scale the coords
    low_d <- scale(low_d)
    attr(low_d, 'scaled:center') <- NULL
    attr(low_d, 'scaled:scale') <- NULL
    err <- forward_cost(low_d, true_dist, std = std)

    return(list(par = low_d, value = err, iters = iter))
}
