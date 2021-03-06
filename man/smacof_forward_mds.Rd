% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forward_methods.R
\name{smacof_forward_mds}
\alias{smacof_forward_mds}
\title{Do forward MDS by using the SMACOF algorithm defined on page 191 of Borg and Groenen.}
\usage{
smacof_forward_mds(high_d, weights, dist.func = euclidean.dist,
  thresh = 1e-05, max.iters = 1000, n.inits = 10, seed = NULL,
  std = TRUE, symm = FALSE)
}
\arguments{
\item{high_d}{The high dimensional data of which a low dimensional representation is desired, an n by p matrix where rows represent observations}

\item{dist.func}{The distance function to be used for both low and high D distance computation.}

\item{thresh}{The threshold below which stress must fall before computation ends.}

\item{max.iters}{The maximum number of SMACOF iterations (i.e., the maximum number of times we solve the quadratic system described in chapter 8 of Borg and Groenen) per initialization.}

\item{n.inits}{The number of times the SMACOF algorithm is run from a random configuration. This can be important, as the cost surface is highly nonconvex.}

\item{seed}{Random seed used for initialization}

\item{symm}{Boolean, is the distance function symmetric? We can save on computation if so.}

\item{low_d}{The low_d solution, an n by 2 matrix, the cost of which is to be evaluated.}
}
\value{
List with $par, the optimal configuration as an n by 2 matrix, and $value as the stress of this configuration
}
\description{
Do forward MDS by using the SMACOF algorithm defined on page 191 of Borg and Groenen.
}
