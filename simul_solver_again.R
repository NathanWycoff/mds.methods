#!/usr/bin/Rscript
#  simul_solver_again.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.09.2017
source('forward_methods.R')
source("inverse_methods.R")

w_step <- function(low_d, last_low_d, dist.func) {
    low_d_dist <- good.dist(low_d, dist.func)
    approx_sol <- function(weights) make.B(high_d_dist, last_low_d) %*% last_low_d

}
