context("test-forward.R")

source("../../forward_methods.R")

#Create some data
seed <- 123
set.seed(seed)
p <- 30
n <- 6
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL

test_that("Dumb MDS works", {
    n.inits <- 10
    dist.func <- euclidean.dist

    #Try the algo with two different seeds and see if they
    # converge to the same thing
    weights <- rep(1,p)
    low_d_result <- forward_mds(high_d, weights, dist.func, n.inits = n.inits, 
                                seed = seed)
    low_d <- low_d_result$par

    par(mfrow=c(2,1))
    plot(low_d, lwd = 0, main = paste('MDS baby', 'stress:',low_d_result$value), 
         lty=0)
    text(low_d[,1], low_d[,2], 1:n)

    old_low_d <- low_d

    low_d_result <- forward_mds(high_d, weights, dist.func, n.inits = n.inits,
                               seed =  seed + 1)
    low_d <- low_d_result$par

    plot(low_d, lwd = 0, main = paste('MDS baby', 'stress:',low_d_result$value),
         lty=0)
    text(low_d[,1], low_d[,2], 1:n)

    #Get the two distance matrices and compare them
    low_1 <- as.matrix(dist(low_d))
    low_2 <- as.matrix(dist(old_low_d))
    diff <- norm(low_1 - low_2)
    expect_lte(diff, 0.5)

})
