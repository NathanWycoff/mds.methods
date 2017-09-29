euclidean.dist <- function(x, y) sqrt(sum((x-y)^2))
inner.prod.dist <- function(x, y) t(x) %*% y

#A distance function that:
#1) Returns a matrix instead of a dist object (I need the diagonals for inner prod distance)
#2) Allows for custom distance function specification.
good.dist <- function(X, dist.func) {
    n <- dim(X)[1]
    dists <- matrix(0, ncol = n, nrow = n)
    for (i in 1:n) {
        for (j in 1:n) {
            dists[i,j] <- dist.func(X[i,], X[j,])
        }
    }
    return(dists)
}

forward_cost <- function(low_d, high_d_dist, k, dist.func) {
    low_d <- matrix(low_d, ncol = k)
    low_d_dist <- good.dist(low_d, dist.func)
    diff_mat <- low_d_dist - high_d_dist

    #Calculate the stress of all elements
    stress <- sum(diff_mat^2)

    return(stress)
}


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
                   thresh = 1e-5, max.iters = 1000, n.inits = 10) {
    #Center high D distances
    high_d <- scale(high_d)
    high_d <- high_d %*% diag(weights)

    #Get our high D distance
    true_dist <- good.dist(high_d, dist.func)

    #Run a bunch of smacof algos
    results <- lapply(1:n.inits, function(i)
                    single_smacof(true_dist, dist.func = dist.func,
                        thresh = thresh, max.iters = max.iters))

    #Get lowest cost result
    costs <- sapply(results, function(i) i$cost)

    return(results[[which.min(costs)]])
}

single_smacof <- function(true_dist, dist.func = euclidean.dist,
                   thresh = 1e-5, max.iters = 1000) {
    #Make a random initial lowD matrix
    low_d <- matrix(rnorm(n*2), ncol = 2)

    ##Do a SMACOF iteration
    diff <- Inf
    last_err <- Inf
    iter <- 0
    while (diff > thresh && iter < max.iters) {
        iter <- iter + 1

        #Create the B matrix
        B <- make.B(true_dist, low_d)

        #Get our next Z
        low_d <- 1/n * B %*% low_d

        #Print current error
        err <- forward_cost(low_d, true_dist, 2, dist.func)

        #Get difference in error
        diff <- abs(last_err - err)
        last_err <- err
    }

    if (iter == max.iters) {
        print("Warning: Did not converge")
    }

    return(list(low_d = low_d, cost = err))
}
