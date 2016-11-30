#require("irlba")
#require("rARPACK")

#' Calculate maximal eigenvalue for big matrices using singular
#' value decomposition
#'
#' @param X Matrix to be evaluated (can be a bigmatrix)
#' @param method One of 'irlba' or 'rARPACK'
#'
#' @return The maximal eigenvalue 
#' 
#' @import irlba
#' @import rARPACK
#'
bigeigen <- function(X, method = "rARPACK") {
  if (method == "irlba") {
    # Doesn't work currently, not sure why
    matmul <- function(A, B) {
      # Bigalgebra requires matrix/vector arguments
      if (is.null(dim(B)))
        B = cbind(B)
      return(cbind( (A %*% B)[]))
    }
    
    max.eigen = (irlba(X, 1, 1, mult = matmul)$d[1])^2
  } else {
    max.eigen = (svds(X, 1)$d[1])^2
  }
  
  return(max.eigen)
}

#' Fused lasso optimisation with the proximal-gradient method (Chen et al. 2010).
#' Internal method
#'
#' @param X matrix of covariates (n by p)
#' @param Y vector of responses (length n)
#' @param groups vector of group indicators (length n)
#' @param group.names group names (length k)
#' @param XX t(X) %*% X
#' @param XY t(X) %*% Y
#' @param p Number of covariates
#' @param k Number of groups
#' @param samp.sizes Sample size for each group
#' @param lambda Sparsity hyperparameter
#' @param gamma Fusion hyperparameter
#' @param C 
#' @param G 
#' @param epsilon 
#' @param tol 
#' @param num.it 
#' @param lam.max 
#' @param c.flag 
#' @param intercept 
#' @param penalty.factors 
#' @param X.list 
#'
#' @return
#' 
#' @useDynLib fuser
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
genFusedLassoProximal <- function(X, Y, groups, group.names=sort(unique(groups)), 
                                  XX=NULL, XY, p, k,
                                  samp.sizes, lambda, gamma, C, G,
                                  epsilon = 1e-04 * p * dim(C)[2],
                                  tol = 1e-06, num.it = 1000,
                                  lam.max = NULL, c.flag = FALSE, 
                                  intercept = TRUE,
                                  penalty.factors = NULL, X.list = NULL) {

  n.c = dim(C)[2]  # Number of constraints, including sparsity
  
  D = 0.5 * p * n.c
  
  mu = epsilon/(2 * D)  # smoothness parameter
  
  C.t = t(C)
  
  if (is.null(lam.max)) {
    # cat('calculate')
    
    if (!is.null(XX)) {
      lam.max = sum(sapply(XX, function(x) max(eigen(x, symmetric = TRUE, only.values = TRUE)$values)))  # maximal eigenvalue
    } else {
      lam.max = sum(sapply(1:length(group.names), function(g) bigeigen(X.list[[g]])))  # maximal eigenvalue
    }
    
    G.temp = G
    G.temp[lower.tri(G)] = 0
    
    L_U = lam.max + (lambda^2 + 2 * gamma^2 * max(colSums(G.temp)))/mu
  } else {
    G.temp = G
    G.temp[lower.tri(G)] = 0
    L_U = lam.max + (lambda^2 + 2 * gamma^2 * max(colSums(G.temp)))/mu
  }
  
  gc()
  
  L_U.inv = 1/L_U
  
  
  if (intercept) {
    W = matrix(0, p + 1, k)
    B.old = matrix(0, p + 1, k)
    weighted.delta.f = matrix(0, p + 1, k)
  } else {
    W = matrix(0, p, k)
    B.old = matrix(0, p, k)
    weighted.delta.f = matrix(0, p, k)
  }
  
  i = 1
  
  for (i in 1:num.it) {
    
    # Allow for scaling of penalty on non-intercept variables (not applying to fusion)
    if (!is.null(penalty.factors)) {
      B.sparsity = B.old[1:p, ] * matrix(penalty.factors, p, k)
    } else {
      B.sparsity = B.old[1:p, ]
    }
    
    # Need to exclude the intercept variable from the B matrix to avoid penalization
    if (intercept) {
      B.sparsity = rbind(B.sparsity, 0)
      B.fusion = rbind(B.old[1:p, ], 0)
    } else {
      B.fusion = B.old
    }
    
    A.star = cbind(B.sparsity %*% C[, 1:k], B.fusion %*% C[, (k + 1):dim(C)[2]])/mu
    
    # A.star = (B.star %*% C) / mu
    A.star[A.star > 1] = 1
    A.star[A.star < -1] = -1
    
    # Derivative is zero when taken with respect to each beta[,k.i]
    if (!is.null(XX)) {
      delta.lik = lapply(1:k, function(k.i) (XX[[k.i]] %*% B.old[, k.i] - XY[[k.i]])/samp.sizes[k.i])
    } else if (c.flag) {
      delta.lik = lapply(1:k, function(k.i) {
        g = group.names[k.i]
        (doubleCrossProd(X.list[[k.i]], B.old[, k.i, drop = FALSE]) - crossprod(X.list[[k.i]], Y[groups == g]))/samp.sizes[k.i]
      })
    } else {
      delta.lik = lapply(1:k, function(k.i) {
        g = group.names[k.i]
        (crossprod(X.list[[k.i]], X.list[[k.i]] %*% B.old[, k.i]) - crossprod(X.list[[k.i]], Y[groups == g]))/samp.sizes[k.i]
      })
    }
    
    delta.f = do.call("cbind", delta.lik) + A.star %*% C.t
    
    B.new = W - L_U.inv * delta.f
    
    weighted.delta.f = weighted.delta.f + L_U.inv * 0.5 * i * delta.f
    
    Z = -weighted.delta.f
    
    W = (i * B.new + 2 * Z)/(i + 2)
    
    i = i + 1
    
    improvement = sum(abs(B.old - B.new))
    
    cat(improvement, '\n')
    
    if (improvement < tol*p)
      break
    
    B.old = B.new
    
    
  }
  
  if (i >= num.it) 
    warning("Reached max iterations without convergence.")
  
  return(B.new)
}


#' Fused lasso optimisation with proximal-gradient method.
#' (Chen et al. 2010)
#'
#' @param X matrix of covariates (n by p)
#' @param Y vector of responses (length n)
#' @param groups vector of group indicators (length n)
#' @param lambda Sparsity hyperparameter (accepts scalar value only)
#' @param gamma Fusion hyperparameter (accepts scalar value only)
#' @param G Matrix of pairwise group information sharing weights (K by K)
#' @param mu Smoothness parameter for proximal optimization
#' @param tol Tolerance for optimization
#' @param num.it Number of iterations
#' @param lam.max Maximal eigenvalue of t(X) %*% X (will be calculate 
#' if not provided)
#' @param c.flag Whether to use Rcpp to calculate cross-products (see Details). 
#' @param intercept Whether to include a (group-specific) intercept term.
#' @param penalty.factors Weights for sparsity penalty.
#' 
#' @details 
#' The proximal algorithm uses t(X)%*%X and t(X)%*%Y. The function will attempt to 
#' pre-calculate these values to speed up computation. This may not always be possible due to
#' memory restrictions; at present this is only done for p < 10,000. When p > 10,000, 
#' crossproducts are calculated explicitly; calculation can be speeded up by using
#' Rcpp code (setting c.flag=TRUE).
#'
#' @return
#' A matrix with the linear coefficients for each group (p by k).
#' 
#' @export
#'
#' @examples
#' set.seed(123)
#' # Generate simple heterogeneous dataset
#' k = 4 # number of groups
#' p = 100 # number of covariates
#' n.group = 15 # number of samples per group
#' sigma = 0.05 # observation noise sd
#' groups = rep(1:k, each=n.group) # group indicators
#' # sparse linear coefficients
#' beta = matrix(0, p, k)
#' nonzero.ind = rbinom(p*k, 1, 0.025/k) # Independent coefficients
#' nonzero.shared = rbinom(p, 1, 0.025) # shared coefficients
#' beta[which(nonzero.ind==1)] = rnorm(sum(nonzero.ind), 1, 0.25) 
#' beta[which(nonzero.shared==1),] = rnorm(sum(nonzero.shared), -1, 0.25)
#' 
#' X = lapply(1:k, function(k.i) matrix(rnorm(n.group*p),n.group, p)) # covariates 
#' y = sapply(1:k, function(k.i) X[[k.i]] %*% beta[,k.i] + rnorm(n.group, 0, sigma)) # response
#' X = do.call('rbind', X)
#' 
#' # Pairwise Fusion strength hyperparameters (tau(k,k'))
#' # Same for all pairs in this example
#' G = matrix(1, k, k) 
#' 
#' # Use L1 fusion to estimate betas (with near-optimal sparsity and 
#' # information sharing among groups)
#' beta.estimate = fusedLassoProximal(X, y, groups, lambda=0.01, tol=3e-3, 
#'                                    gamma=0.01, G, intercept=FALSE,
#'                                    num.it=500) 
fusedLassoProximal <- function(X, Y, groups, lambda, gamma, G, mu = 1e-04, tol = 1e-06, num.it = 1000, 
                               lam.max = NULL, c.flag = FALSE, intercept = TRUE, 
                               penalty.factors=NULL) {
  
  group.names = sort(unique(groups))
  samp.sizes = table(groups)[group.names]
  
  k = length(group.names)
  p = dim(X)[2]
  
  if(is.null(penalty.factors)) 
    penalty.factors = rep(1,dim(X)[2])
  
  if (intercept) {
    X = cbind(X, 1)
  }
  
  # Rough limit for keeping everything in memory
  if (p < 10000) {
    XX.list = lapply(group.names, function(x) {
      indices = groups == x
      crossprod(X[indices, ])
    })
    
    XY.list = lapply(group.names, function(x) {
      indices = groups == x
      crossprod(X[indices, ], Y[indices])
    })
    X.list = NULL
  } else {
    # Otherwise calculate on the fly (takes longer, but saves memory)
    XX.list = NULL
    XY.list = NULL
    
    # Avoid indexing over large matrix
    X.list = list()
    
    for (k.i in 1:k) {
      X.list[[k.i]] = X[groups == group.names[k.i], ]
    }
     
    X = NULL # No longer need full matrix, delete to save memory
  }
  
  
  C.fusion = matrix(0, k, choose(k, 2))
  
  edge.i = 1
  
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      C.fusion[i, edge.i] = G[i, j]
      C.fusion[j, edge.i] = -G[i, j]
      edge.i = edge.i + 1
    }
  }
  
  C.sparsity = lambda * diag(k)
  
  C = cbind(C.sparsity, gamma * C.fusion)
  # cat(dim(C), '\n')
  
  return(genFusedLassoProximal(X, Y, groups, group.names, XX.list, XY.list, p, k, samp.sizes, lambda, gamma, C, G, epsilon = mu * p * dim(C)[2], tol, num.it, lam.max, 
                               c.flag=c.flag, intercept = intercept, penalty.factors = penalty.factors, X.list = X.list))
}
