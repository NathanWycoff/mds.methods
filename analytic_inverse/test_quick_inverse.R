source('quick_inverse.R')
source('forward_methods.R')

##Generate some weights and try to recover them.
#Some params
n <- 20#How many points?
p <- 6#How many high D dims?
k <- 2#How many low D dims (never not 2, but I like to be overly general)
n.inits <- 1
dist.func <- function(x, y, w = 1) sqrt(sum(w*(x - y)^2))
seed <- 123


#Generate the weights and induced lowD points
high_d <- matrix(rnorm(n*p), ncol = p)
user_weights <- rgamma(p, 1, 1)
user_weights <- user_weights / sum(user_weights)
low_d_result <- forward_mds(high_d, k, user_weights, dist.func, n.inits, seed + 1)
user_low_d <- low_d_result$par

#Some inital high d points
#init_high_d <- matrix(rnorm(n*p), ncol = p)
init_high_d <- high_d

#Precompute some distance matrices
DELTA <- as.matrix(dist(user_low_d))
D <- as.matrix(dist(init_high_d))
Z <- init_high_d
U <- user_low_d

#Try to recover those weights based on the low D points
weights <- user_weights
learn_rate <- 0.01
diff <- as.numeric('Inf')
thresh <- 1e-5

while (diff > thresh) {
    grad <- maj_grad(weights, U, DELTA, D, Z, dist.func)
    old_weights <- weights
    weights <- weights - learn_rate * grad
    weights[weights < 0] <- 0
    diff <- norm(old_weights - weights, '2')
    print(weights)
}
