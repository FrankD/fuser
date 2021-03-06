% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l1_fusion.R
\name{fusedLassoProximal}
\alias{fusedLassoProximal}
\title{Fused lasso optimisation with proximal-gradient method.
(Chen et al. 2010)}
\usage{
fusedLassoProximal(X, Y, groups, lambda, gamma, G, mu = 1e-04, tol = 1e-06,
  num.it = 1000, lam.max = NULL, c.flag = FALSE, intercept = TRUE,
  penalty.factors = NULL, conserve.memory = p >= 10000, scaling = TRUE)
}
\arguments{
\item{X}{matrix of covariates (n by p)}

\item{Y}{vector of responses (length n)}

\item{groups}{vector of group indicators (length n)}

\item{lambda}{Sparsity hyperparameter (accepts scalar value only)}

\item{gamma}{Fusion hyperparameter (accepts scalar value only)}

\item{G}{Matrix of pairwise group information sharing weights (K by K)}

\item{mu}{Smoothness parameter for proximal optimization}

\item{tol}{Tolerance for optimization}

\item{num.it}{Number of iterations}

\item{lam.max}{Maximal eigenvalue of \code{t(X) \%*\% X} (will be calculate
if not provided)}

\item{c.flag}{Whether to use Rcpp for certain calculations (see Details).}

\item{intercept}{Whether to include a (group-specific) intercept term.}

\item{penalty.factors}{Weights for sparsity penalty.}

\item{conserve.memory}{Whether to calculate XX and XY on the fly, conserving memory
at the cost of speed. (True by default iff p >= 10000)}

\item{scaling}{if TRUE, scale the sum-squared loss for each group by 1/n_k
where n_k is the number of samples in group k.}
}
\value{
A matrix with the linear coefficients for each group (p by k).
}
\description{
Fused lasso optimisation with proximal-gradient method.
(Chen et al. 2010)
}
\details{
The proximal algorithm uses \code{t(X) \%*\% X} and \code{t(X) \%*\% Y}. The function will attempt to
pre-calculate these values to speed up computation. This may not always be possible due to
memory restrictions; at present this is only done for p < 10,000. When p > 10,000,
crossproducts are calculated explicitly; calculation can be speeded up by using
Rcpp code (setting c.flag=TRUE).
}
\examples{
set.seed(123)
# Generate simple heterogeneous dataset
k = 4 # number of groups
p = 100 # number of covariates
n.group = 15 # number of samples per group
sigma = 0.05 # observation noise sd
groups = rep(1:k, each=n.group) # group indicators
# sparse linear coefficients
beta = matrix(0, p, k)
nonzero.ind = rbinom(p*k, 1, 0.025/k) # Independent coefficients
nonzero.shared = rbinom(p, 1, 0.025) # shared coefficients
beta[which(nonzero.ind==1)] = rnorm(sum(nonzero.ind), 1, 0.25)
beta[which(nonzero.shared==1),] = rnorm(sum(nonzero.shared), -1, 0.25)

X = lapply(1:k,
           function(k.i) matrix(rnorm(n.group*p),
                                n.group, p)) # covariates
y = sapply(1:k,
           function(k.i) X[[k.i]] \%*\% beta[,k.i] +
                           rnorm(n.group, 0, sigma)) # response
X = do.call('rbind', X)

# Pairwise Fusion strength hyperparameters (tau(k,k'))
# Same for all pairs in this example
G = matrix(1, k, k)

# Use L1 fusion to estimate betas (with near-optimal sparsity and
# information sharing among groups)
beta.estimate = fusedLassoProximal(X, y, groups, lambda=0.01, tol=3e-3,
                                   gamma=0.01, G, intercept=FALSE,
                                   num.it=500)
}
