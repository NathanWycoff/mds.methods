#!/usr/bin/Rscript
#  inverse_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.03.2017
source('forward_methods.R')
source("inverse_methods.R")
source('simul_solver.R')

##Generate some weights and try to recover them.
#Some params
thresh <- 1e-5
max.iters <- 1000

#Create some data
seed <- 12345678
set.seed(seed)
p <- 5
n <- 5
high_d <- matrix(rnorm(n*p), ncol = p)
high_d <- scale(high_d)
attr(high_d, 'scaled:center') <- NULL
attr(high_d, 'scaled:scale') <- NULL
k <- 2#How many low D dims (never not 2, but I like to be overly general)

#Specify some params
n.inits <- 10
dist.func <- euclidean.dist


#####Try to get a slight perturbation of the weights
#Generate the 'initial' data
set.seed(seed)
init_weights <- rgamma(p, 1, 1)
init_weights <- init_weights / sum(init_weights)
low_d_result <- forward_mds(high_d, k, init_weights, dist.func, n.inits, seed + 1)
init_low_d <- low_d_result$par

#Different weights
#perturb_weights <- rep(1,p)
#perturb_weights[2] <- 10
set.seed(seed+1)
perturb_weights <- rgamma(p,1,1)
perturb_weights <- perturb_weights / sum(perturb_weights)
low_d_result <- forward_mds(high_d, k, perturb_weights, dist.func, n.inits, seed + 1)
perturb_low_d <- low_d_result$par

#Plot the two results
par(mfrow=c(2,1))
plot(init_low_d, main='init', lty=0, lwd=0)
text(init_low_d[,1], init_low_d[,2], 1:n)
plot(perturb_low_d, main='perturb', lty=0, lwd=0)
text(perturb_low_d[,1], perturb_low_d[,2], 1:n)

###############################################
###Confirm that SMACOF approx and exact cost are very close if there hasnt been perturbation
#low_d_dist <- good.dist(init_low_d, euclidean.dist)
#a <- approx_smacof_inverse_cost(init_weights, low_d_dist, 
#                                 high_d, init_low_d, n.inits, dist.func)
#
#b <- smacof_inverse_cost(init_weights, low_d_dist, high_d, k, n.inits, dist.func) 
#
#print(paste("Approx SMACOF cost with no perturb:", a))
#print(paste("Exact SMACOF cost with no perturb:", b))
#
################################################
##If we start with exactly the correct weights, how far do we move?
#forward.n.inits <- 10
#system.time(optim_result <- approx_smacof_inverse_step(init_low_d, init_low_d, 
#                                          high_d, init_weights, dist.func, 
#                                          forward.n.inits, seed + 1))
#infer_weights <- optim_result$par / sum(optim_result$par)
#
#print("Init Weights:")
#print(init_weights)
#print("Perturb Weights:")
#print(perturb_weights)
#print("My guess:")
#print(infer_weights)
#print(paste('mine',sum(abs(infer_weights - perturb_weights))))
#print(paste('init',sum(abs(init_weights - perturb_weights))))
#print(paste('unif',sum(abs(mean(perturb_weights) - perturb_weights))))


################################################
###If we start with the correct weights and only perturb a little bit, can we recover the new weights using only 1 inverse step?
#forward.n.inits <- 10
#system.time(optim_result <- approx_smacof_inverse_step(perturb_low_d, init_low_d, 
#                                          high_d, init_weights, dist.func, 
#                                          forward.n.inits, seed + 1))
#infer_weights <- optim_result$par / sum(optim_result$par)
#
#print("Init Weights:")
#print(init_weights)
#print("Perturb Weights:")
#print(perturb_weights)
#print("My guess:")
#print(infer_weights)
#print(paste('mine',sum(abs(infer_weights - perturb_weights))))
#print(paste('init',sum(abs(init_weights - perturb_weights))))
#print(paste('unif',sum(abs(mean(perturb_weights) - perturb_weights))))

#################################################
## Run the entire inverse algo, trying to recover the weights
system.time(result <- approx_smacof_simul(perturb_low_d, init_low_d, high_d, 
                               init_weights, dist.func, forward.n.inits, thresh, 
                               max.iters-1, seed))
infer_weights <- result$weights

print("Init Weights:")
print(init_weights)

print("Perturb Weights:")
print(perturb_weights)

print("My guess:")
print(infer_weights)

print(paste('mine',sum(abs(infer_weights - perturb_weights))))
print(paste('init',sum(abs(init_weights - perturb_weights))))
print(paste('unif',sum(abs(mean(perturb_weights) - perturb_weights))))


##Get the SMACOF cost of the two solutions
#p_dist <- good.dist(perturb_low_d, dist.func)
#approx_smacof_inverse_cost(infer_weights, p_dist, 
#                                 high_d, init_low_d, n.inits, dist.func)
#approx_smacof_inverse_cost(init_weights, p_dist, 
#                                 high_d, init_low_d, n.inits, dist.func)
#approx_smacof_inverse_cost(perturb_weights, p_dist, 
#                                 high_d, init_low_d, n.inits, dist.func)

par(mfrow=c(2,2))
plot(result$Zs[[999]], main='999', lty=0, lwd=0)
text(result$Zs[[999]][,1],result$Zs[[999]][,2],  1:n)
plot(result$Zs[[998]], main='998', lty=0, lwd=0)
text(result$Zs[[998]][,1],result$Zs[[998]][,2],  1:n)
plot(perturb_low_d, main='perturb', lty=0, lwd=0)
text(perturb_low_d[,1], perturb_low_d[,2], 1:n)
plot(init_low_d, main='init', lty=0, lwd=0)
text(init_low_d[,1], init_low_d[,2], 1:n)
