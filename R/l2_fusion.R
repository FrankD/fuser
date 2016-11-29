# L2 Fusion approach optimization
#library("glmnet")
#library("Matrix")

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
#' @examples
#' add(1, 1)
#' add(10, 1)
#' 
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
#' @export
#'
#' @examples
fusedL2DescentGLMNet <- function(new.x, new.x.f, new.y, groups, lambda, gamma=1,
                                 ...) {
  
  # Incorporate fusion penalty global hyperparameter
  new.x.f = new.x.f * sqrt(gamma * 
                             (dim(new.x)[1] + dim(new.x.f)[1]))
  new.x = rBind(new.x, new.x.f)
  
  group.names = sort(unique(groups))
  num.groups = length(group.names)
  
  glmnet.result = glmnet(new.x, new.y, standardize=FALSE, intercept=FALSE, ...)
  
  beta.mat = array(NA, c(dim(new.x)[2]/num.groups, num.groups, length(lambda)))
  
  for(lambda.i in 1:length(lambda)) {
    
    coef.temp = coef(glmnet.result, 
                     s=lambda[lambda.i]*length(groups)/dim(new.x)[1]) # Correction for extra dimensions
    beta.mat[,,lambda.i] = coef.temp[2:length(coef.temp)]
  }
  
  return(beta.mat)
}
