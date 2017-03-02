# L2 Fusion approach optimization


# 

#' Generate block diagonal matrices to allow for fused L2 optimization with glmnet.
#' 
#' @param X covariates matrix (n by p).
#' @param response vector (length n).
#' @param groups vector of group indicators (ideally factors, length n)
#' @param G matrix representing the fusion strengths between pairs of 
#' groups (K by K). Zero entries are assumed to be independent pairs.
#' @param intercept whether to include an (per-group) intercept in the model
#' @param penalty.factors vector of weights for the penalization of 
#' each covariate (length p)
#' 
#' @return A list with components X, Y, X.fused and penalty, where
#' X is a n by pK block-diagonal bigmatrix, Y is a 
#' re-arranged bigvector of length n, and X.fused is a 
#' choose(K,2)*p by pK bigmatrix encoding the fusion penalties.
#' 
#' @import Matrix
#'  
#' @examples
#' set.seed(123)
#' 
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
#' # Generate block diagonal matrices
#' transformed.data = generateBlockDiagonalMatrices(X, y, groups, G)
generateBlockDiagonalMatrices <- function(X, Y, groups, G, intercept=FALSE,
                                          penalty.factors=rep(1, dim(X)[2])) {
  group.names = sort(unique(groups))
  num.groups = length(group.names)
  
  num.pairs = num.groups*(num.groups-1)/2
  
  if(intercept) X = cbind(X, matrix(1, dim(X)[1], 1)) # Include intercept
  
  new.y = Matrix(0, length(Y)+num.pairs*dim(X)[2],1, sparse=TRUE)
  new.x = Matrix(0, dim(X)[1], dim(X)[2]*num.groups,
                 sparse=TRUE)
  new.x.f = Matrix(0, num.pairs*dim(X)[2], dim(X)[2]*num.groups,
                   sparse=TRUE)
  
  row.start = 1
  col.start = 1
  
  # Add data into block matrix
  for(group.i in 1:num.groups) {
    group.inds = groups==group.names[group.i]
    
    row.range = row.start:(row.start+sum(group.inds)-1)
    col.range = col.start:(col.start+dim(X)[2]-1)
    
    new.y[row.range] = Y[group.inds]
    new.x[row.range, col.range] = X[group.inds,]
    
    row.start = row.start + sum(group.inds)
    col.start = col.start + dim(X)[2]
  }
  
  col.start.i = 1
  row.start = 1
  
  # Add gamma contraints into block matrix
  for(group.i in 1:(num.groups-1)) {
    
    col.start.j = col.start.i + dim(X)[2]
    
    col.range.i = col.start.i:(col.start.i+dim(X)[2]-1)
    
    for(group.j in (group.i+1):num.groups) {
      
      tau = G[group.i, group.j]
      
      row.range = row.start:(row.start+dim(X)[2]-1)
      col.range.j = col.start.j:(col.start.j+dim(X)[2]-1)
      
      new.x.f[cbind(row.range, col.range.i)] = sqrt(tau)
      new.x.f[cbind(row.range, col.range.j)] = -sqrt(tau)
      
      if(intercept) { # Don't fuse intercept
        new.x.f[row.range[length(row.range)], 
                col.range.i[length(col.range.i)]] = 0
        new.x.f[row.range[length(row.range)], 
                col.range.j[length(col.range.j)]] = 0
      }
      
      row.start = row.start + dim(X)[2]
      col.start.j = col.start.j + dim(X)[2]
    }
    
    col.start.i = col.start.i + dim(X)[2]
  }
  
  if(intercept) {
    penalty = rep(c(penalty.factors, 0), num.groups)
  } else {
    penalty = rep(penalty.factors, num.groups)
  }
  
  return(list(X=new.x, Y=new.y, X.fused=new.x.f, penalty=penalty))
}

#' Optimise the fused L2 model with glmnet (using transformed input data) 
#'
#' @param transformed.x Transformed covariates (output of generateBlockDiagonalMatrices)
#' @param transformed.x.f Transformed fusion constraints (output of generateBlockDiagonalMatrices)
#' @param transformed.y Transformed response (output of generateBlockDiagonalMatrices)
#' @param num.groups Number of groups (K)
#' @param group.names Names for the groups (optional)
#' @param lambda.sparse Sparsity penalty hyperparameter
#' @param ... Further options passed to glmnet.
#'
#' @return Matrix of fitted beta values.
#' @import glmnet
#' @export
#'
#' @return 
#' A matrix with the linear coefficients for each group (p by k).
#'
#' @examples
#' 
#' #' set.seed(123)
#' 
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
#' # Generate block diagonal matrices
#' transformed.data = generateBlockDiagonalMatrices(X, y, groups, G)
#' 
#' # Use L2 fusion to estimate betas (with near-optimal information sharing among groups)
#' beta.estimate = fusedL2DescentGLMNet(transformed.data$X, transformed.data$X.fused,
#'                                      transformed.data$Y, groups, lambda=c(0,0.001,0.1,1),
#'                                      gamma=0.001)
fusedL2DescentGLMNet <- function(new.x, new.x.f, new.y, groups, lambda, gamma=1,
                                 ...) {
  
  # Incorporate fusion penalty global hyperparameter
  new.x.f = new.x.f * sqrt(gamma * 
                             (dim(new.x)[1] + dim(new.x.f)[1]))
  new.x = rBind(new.x, new.x.f)
  
  group.names = sort(unique(groups))
  num.groups = length(group.names)
  
  glmnet.result = glmnet(new.x, new.y, standardize=FALSE, ...)
  
  beta.mat = array(NA, c(dim(new.x)[2]/num.groups, num.groups, length(lambda)))
  
  for(lambda.i in 1:length(lambda)) {
    
    coef.temp = coef(glmnet.result, 
                     s=lambda[lambda.i]*length(groups)/dim(new.x)[1]) # Correction for extra dimensions
    beta.mat[,,lambda.i] = coef.temp[2:length(coef.temp)]
  }
  
  return(beta.mat)
}

# Improvement for the general fused L2 formulation
likelihoodImprovementGFL2R <- function(beta.new, beta.old, X, Y, groups, group.names,
                                       lambda.sparse, lambda.fused, G) {
  
  fusion.penalty.old = 0
  fusion.penalty.new = 0
  
  for(i in 1:(dim(beta.new)[2]-1)) {
    for(j in (i+1):dim(beta.new)[2]) {
      fusion.penalty.new = fusion.penalty.new + G[i,j] * sum((beta.new[,i]-beta.new[,j])^2)
      fusion.penalty.old = fusion.penalty.old + G[i,j] * sum((beta.old[,i]-beta.old[,j])^2)
    }
  }
  
  lh.new = 0
  lh.old = 0
  
  for(group.i in 1:dim(beta.new)[2]) {
    group.name = group.names[group.i]
    
    lh.new = lh.new + sum((Y[groups==group.name] - X[groups==group.name,] %*% 
                             beta.new[,group.i])^2)
    lh.old = lh.old + sum((Y[groups==group.name] - X[groups==group.name,] %*% 
                             beta.old[,group.i])^2)
  }
  
  lh.new = 0.5*lh.new + lambda.sparse*sum(abs(beta.new)) + 
    lambda.fused*fusion.penalty.new
  lh.old = 0.5*lh.old + lambda.sparse*sum(abs(beta.old)) + 
    lambda.fused*fusion.penalty.old
  
  #if(lh.new > lh.old) stop('Likelihood did not decrease.')
  
  return(list(improvement=lh.old - lh.new, lh.old=lh.old, lh.new=lh.new))
}

# Update one row of the beta matrix (i.e. the inseparable part of the greater 
# separable problem) with L2 fused. Following Friedman et al (2007)
betaRowFusedL2Update <- function(row.i, X, Y, V.inv, XVX, XVY, gamma,
                                 G, beta.mat, groups, group.names) {
  K = dim(beta.mat)[2] # Number of groups
  
  gamma_2 = 2*gamma
  
  for(beta.i in 1:K) {
    beta.row = beta.mat[row.i,]
    
    # Least squares estimate without additional penalty beyond L2 fusion
    if(!is.null(XVX)) {
      ls.beta = (XVY[[beta.i]][row.i] - XVX[[beta.i]][row.i,-row.i] %*% beta.mat[-row.i,beta.i] +
                   gamma_2*G[beta.i,-beta.i]%*%beta.row[-beta.i]) / 
        (XVX[[beta.i]][row.i, row.i] + gamma_2*sum(G[beta.i,-beta.i]))
    } else {
      group.inds  = groups == group.names[beta.i]
      ls.beta = (crossprod(X[group.inds,row.i], Y[group.inds]) - 
                   crossprod(X[group.inds,row.i], X[group.inds,-row.i] %*% beta.mat[-row.i,beta.i]) +
                   gamma_2*G[beta.i,-beta.i]%*%beta.row[-beta.i]) / 
        (crossprod(X[group.inds, row.i]) + gamma_2*sum(G[beta.i,-beta.i]))
    }
    
    # Now soft-threshold to get sparse lasso solution
    #lasso.beta = sign(ls.beta) * max(abs(ls.beta) - lambda.sparse, 0)
    
    beta.mat[row.i,beta.i] = ls.beta
  }
  
  return(beta.mat[row.i,])
}

# Fused L2 penalty with block coordinate descent method (Friedman et al. 2007)
fusedL2DescentDirect <- function(X, Y, V.inv, groups, gamma, G,
                                 beta.init=matrix(0, dim(X)[2], length(unique(groups))),
                                 eps = 1e-6, c.flag=FALSE) {
  p = dim(X)[2] # Dimensionality
  group.names = sort(unique(groups))
  
  if(p < 10000) {
    XVY = lapply(group.names, function(g) crossprod(X[groups==g,], V.inv[[which(group.names==g)]] %*% Y[groups==g]))
    XVX = lapply(group.names, function(g) crossprod(X[groups==g,], V.inv[[which(group.names==g)]] %*% X[groups==g,]))
  } else {
    XVY = NULL
    XVX = NULL
  }
  
  beta.mat = beta.init
  
  if(!c.flag) {
    
    diff = Inf # maximum likelihood update
    
    beta.new = beta.mat
    
    while(diff > eps) {
      
      for(row.i in 1:p) {
        beta.new[row.i,] = betaRowFusedL2Update(row.i, X, Y, V.inv, XVX, XVY, 
                                                gamma, G, beta.new, groups, group.names)
        #cat(beta.new[row.i,], '\n')
      }
      
      # Not the correct likelihood for LMM
      #diff = likelihoodImprovementGFL2R(beta.new, beta.mat, X, Y, groups, group.names,
      #                                  lambda.sparse, lambda.fused, G) 
      
      diff = sqrt(mean((beta.mat-beta.new)^2))
      
      beta.mat = beta.new
    }
  } else {
    # Not implemented
    # Everything done inside C++ routine
    # 
    # if(!is.null(XX)) {
    #   beta.mat = betaFusedL2Optim(X, Y, groups, XX, XY, lambda.sparse, 
    #                               lambda.fused, G, beta.mat, eps)
    # } else {
    #   beta.mat = betaFusedL2Optim2(X, Y, groups, lambda.sparse, 
    #                                lambda.fused, G, beta.mat, eps)
    # }
  }
  
  
  return(beta.mat)
}
