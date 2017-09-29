#################################################################################
##### Functions for forwards
#An inner product function
inner.prod.dist <- function(x, y) t(x) %*% y
euclidean.dist <- function(x, y) sqrt(sum((x-y)^2))

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

forward_mds <- function(high_d, k, weights, dist.func, n.inits, seed) {
    #Cholesky decomp of the diagonal matrix
    weights <- sqrt(weights)
    high_d <- high_d %*% diag(weights)

    #Create the distance matrix
    n <- dim(high_d)[1]
    high_d_dist <- good.dist(high_d, dist.func)
    
    #Randomly Init Points
    set.seed(seed)
    inits <- lapply(1:n.inits, function(x) rnorm(n*k))

    #Get optimal points
    z_optimals <- lapply(inits, function(init) optim(init, forward_cost, method = "BFGS", high_d_dist = high_d_dist, k = k, dist.func = dist.func, control = list('abstol' = 1e-15)))

    optimal <- which.min(lapply(z_optimals, function(z) z$value))
    z_optimal <- z_optimals[[optimal]]

    #Shape it into a matrix
    z_optimal$par <- matrix(z_optimal$par, ncol = k)

    return(z_optimal)
}
