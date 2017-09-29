#This gives the gradient of the majorizing function
#This gradient (MIGHT) bound the gradient of stress at the point the functions are equal.
#@param weights - The weights about which to evaluate the gradient
#@param U - The user specified low D coords
#@param DELTA - low D distances
#@param D - high D distances
#@param Z - The high D points induced by the last best weights
#@param dist.func - distance function of two variates, parameterized by a weight vector

maj_grad <- function(weights, U, DELTA, D, Z, dist.func) {
    n <- dim(U)[1]
    d <- dim(U)[2]
    p <- dim(Z)[2]


    grads <- c()
    for (k in 1:p) {
        grad <- 0

        for (i in 1:n) {
            for (j in 1:(i-1)) {
                for (l in 1:d) {
                    g <- -U[i,l] + Z[j,l]
                    insum <- 0
                    for (a in (1:n)[-j]) {
                        insum <- insum + DELTA[j,a] * (Z[a,l] - Z[j,l])
                        insum <- insum * 1 / D[j,a]^2

                        #Calculate difference in dist using finite diffs so that we can accomadate any distance func
                        h <- rep(0, p)
                        h[k] <- 1e-8
                        dddw <- (dist.func(Z[j, ], Z[a,], weights + h) - dist.func(Z[j, ], Z[a,], weights)) / h[k]

                        insum <- insum * dddw
                    }
                }
                grad <- grad + insum
            }
        }
        grads <- c(grads, grad)
    }

    return(grads)
}
