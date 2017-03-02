#' Fit a linear mixed model (using the nlme package) with
#' fusion over the fixed effects. Only l2 fusion is currently implemented.
#' 
#' @description A linear mixed effects model is fitted for each group using 
#' the lme function of package nlme. Subsequently, the l2 fusion approach is 
#' applied to share information about the fixed effects among the groups.
#' 
#' Unlike the functions for fusion estimation in fixed effect linear models, 
#' this function does not currently scale to high-dimensional data with a 
#' large number of fixed effects. Consequently, l1 penalisation of the 
#' fixed effects has not been implemented because it is assumed that n > p.
#'
#' @param formula A two-sided linear formula object describing the fixed effects,
#' identical to the input of the 'lme' function. Note that we have assumed a single 
#' response variable.
#' @param random A one-sided formula or list of one-sided formulas describing the 
#' random effects. Identical to the input of the 'lme' function.
#' @param data A data frame containing the variables names in fixed and random.
#' @param groups A vector of group indicators (length n)
#' @param gamma The fusion hyperparameter (must be a scalar)
#' @param G Matrix of pairwise group information sharing weights (K by K). 
#' Default is equal sharing across all pairs.
#' 
#' @return A list of class fused.lmm, with elements 'beta', the fixed effects 
#' estimated with fusion penalty information sharing across groups, and 
#' 'lme.fits', the fitted lme models for each group. Note that the lme fits 
#' will also contain estimates of the fixed effects, but estimated without 
#' information sharing.
#' 
#' @import nlme
#' @import mgcv
fusionLMM <- function(formula, random, data, groups, gamma, 
                      G=matrix(1,length(unique(groups)),
                               length(unique(groups)))) {
  
  group.ids = sort(unique(groups))
  
  # Fit lmms
  lme.models = lapply(group.ids, function(group.id) {
    model = lme(formula, random=random, data=data[groups==group.id,],
                na.action=na.omit, control=lmeControl(returnObject=TRUE))
    params = extract.lme.cov(model, data=data[groups==group.id,])
    model$V = params$V; model$G=params$G; model$Z=params$Z
    return(model)
  })
  names(lme.models) = group.ids
  
  V.mats = lapply(lme.models, function(model) model$V)
  
  # get response variable
  response.var = labels(terms(getResponseFormula(formula)))

  X = model.matrix(formula, data)
  Y = data[[response.var]]
  
  # Fit fused l2 model to fixed effects
  beta = refitLMMFixedEffects(X, Y, V.mats, groups, gamma=gamma, G=G)
  colnames(beta) = group.ids
  
  # Update models with beta estimates under fusion
  lme.models = lapply(group.ids, function(group.id) {
    model = lme.models[[group.id]]
    model$beta.fused = beta[,group.id,drop=FALSE]
    model$resid.fixed = Y[groups==group.id] - X[groups==group.id,] %*% model$beta.fused
    return(model)
  })
  names(lme.models) = group.ids
  
  class(lme.models) = 'fused.lmm'
  
  return(lme.models)
}



#' Apply fusion to fixed effects of a linear mixed model. 
#' 
#' @description In order to use 
#' this function, we first require standard fitting of a linear mixed model 
#' for each group and then estimation of the variance-covariance matrices
#' of the response.
#'
#' @param X Fixed effect covariates (rhs of the formula)
#' @param Y Response variable
#' @param V variance-covariance matrices of the response, a list of K n_k 
#' by n_k matrices, where K is the number of groups.
#' @param groups A vector of group indicators (length n)
#' @param gamma The fusion hyperparameter (must be a scalar)
#' @param G Matrix of pairwise group information sharing weights (K by K)
#' @params eps Convergence threshold for the gradient descent algorithm
#' 
#' @return A matrix of fixed effects for each group.
refitLMMFixedEffects <- function(X, Y, V, groups, gamma=0, G, eps=1e-5) {
  
  V.inv = lapply(V, solve)
  
  return(fusedL2DescentDirect(X, Y, V.inv, groups, gamma, G, eps=eps))
  
}

#' Prediction method for fused.lmm objects
#' 
#' @description BLUP estimation of the response for new covariate data.
#'
#' @param object An S3 object of class 'fused.lmm'
#' @param newdata.fixed Model matrix for the fixed effects
#' @param newdata.random Model matrix for the random effects
#' @param groups A vector of group indicators (length n). Note that group 
#' indicators must be identical to the ones used to train the fused.lmm
#' model.
#' 
#' @return A matrix of fixed effects for each group.
predict.fused.lmm <- function(object, newdata.fixed, newdata.random, groups, ...) {
  group.ids = sort(unique(groups))
  
  # Predict for each group
  y.predict = lapply(group.ids, function(group.id) {
    model = object[[group.id]]
    y.predict = newdata.fixed %*% model$beta.fused + newdata.random %*% model$G %*% t(model$Z) %*% 
      solve(model$X) %*% model$resid.fixed
  })
  browser()
  # Something off with calculation of Z -- needs further work
  return(do.call('rbind', y.predict))
  
}

#' Function for extracting Z, G, and V from an lme object.
#' 
#' @description This function has been adapted from the same function in 
#' the package mgcv, to allow for extraction of more parameters. It calculates
#' the variance-covariance matrix Z of the data, and also returns the 
#' design matrix of random effects Z and the covariance matrix of the random 
#' effects G.
#'
#' @param b A fitted object of class 'lme'.
#' @param data The data that we want to calculate Z and V.
#' 
#' @return A list with elements V, Z and G.
#' @import mgcv
extract.lme.cov <- function(b,data) { 
  # function to extract the response data covariance matrix from an lme fitted
  # model object b, fitted to the data in data. "inner" == "finest" grouping 
  # start.level is the r.e. grouping level at which to start the construction, 
  # levels outer to this will not be included in the calculation - this is useful
  # for gamm calculations
  
  start.level=1
  
  if (!inherits(b,"lme")) stop("object does not appear to be of class lme")
  grps<-nlme::getGroups(b) # labels of the innermost groupings - in data frame order
  n<-length(grps)    # number of data
  if (is.null(b$modelStruct$varStruct)) w<-rep(b$sigma,n) ### 
  else 
  { w<-1/nlme::varWeights(b$modelStruct$varStruct) 
  # w is not in data.frame order - it's in inner grouping level order
  group.name<-names(b$groups) # b$groups[[i]] doesn't always retain factor ordering
  order.txt <- paste("ind<-order(data[[\"",group.name[1],"\"]]",sep="")
  if (length(b$groups)>1) for (i in 2:length(b$groups)) 
    order.txt <- paste(order.txt,",data[[\"",group.name[i],"\"]]",sep="")
  order.txt <- paste(order.txt,")")
  eval(parse(text=order.txt))
  w[ind] <- w # into data frame order
  w<-w*b$sigma
  }
  if (is.null(b$modelStruct$corStruct)) V<-diag(n) #*b$sigma^2
  else
  { c.m<-nlme::corMatrix(b$modelStruct$corStruct) # correlation matrices for each group
  if (!is.list(c.m)) V<-c.m
  else
  { V<-matrix(0,n,n)   # data cor matrix
  gr.name <- names(c.m) # the names of the groups
  n.g<-length(c.m)   # number of innermost groups
  j0<-1
  ind<-ii<-1:n
  for (i in 1:n.g) 
  { j1<-j0+nrow(c.m[[i]])-1
  V[j0:j1,j0:j1]<-c.m[[i]]
  ind[j0:j1]<-ii[grps==gr.name[i]]
  j0<-j1+1  
  }
  V[ind,]<-V;V[,ind]<-V # pasting correlations into right place in overall matrix
  # V<-V*b$sigma^2
  }
  }  
  V <- as.vector(w)*t(as.vector(w)*V) # diag(w)%*%V%*%diag(w)
  # ... covariance matrix according to fitted correlation structure
  X<-list()
  grp.dims<-b$dims$ncol # number of Zt columns for each grouping level (inner levels first)
  # inner levels are first in Zt
  Zt<-model.matrix(b$modelStruct$reStruct,data)  # a sort of proto - Z matrix
  # b$groups and cov (defined below have the inner levels last)
  cov<-as.matrix(b$modelStruct$reStruct) # list of estimated covariance matrices (inner level last)
  i.col<-1
  n.levels<-length(b$groups)
  Z<-matrix(0,n,0) # Z matrix
  if (start.level<=n.levels)
  { for (i in 1:(n.levels-start.level+1)) # work through the r.e. groupings inner to outer
  { # get matrix with columns that are indicator variables for ith set of groups...
    # groups has outer levels first 
    if(length(levels(b$groups[[n.levels-i+1]]))==1) { ## model.matrix needs >1 level 
      X[[1]] <- matrix(rep(1,nrow(b$groups))) } else { 
        X[[1]] <- model.matrix(~b$groups[[n.levels-i+1]]-1,
                               contrasts.arg=c("contr.treatment","contr.treatment")) }
    # Get `model matrix' columns relevant to current grouping level...
    X[[2]] <- Zt[,i.col:(i.col+grp.dims[i]-1),drop=FALSE]
    i.col <- i.col+grp.dims[i]
    # tensor product the X[[1]] and X[[2]] rows...
    Z <- cbind(Z,tensor.prod.model.matrix(X))
  } # so Z assembled from inner to outer levels
    # Now construct overall ranef covariance matrix
    Vr <- matrix(0,ncol(Z),ncol(Z))
    start <- 1
    for (i in 1:(n.levels-start.level+1))
    { k <- n.levels-i+1
    for (j in 1:b$dims$ngrps[i]) 
    { stop <- start+ncol(cov[[k]])-1
    Vr[start:stop,start:stop]<-cov[[k]]
    start <- stop+1
    }
    }
    Vr <- Vr*b$sigma^2
    V <- V+Z%*%Vr%*%t(Z)
  }
  
  return(list(V=V, G=Vr, Z=Z))
}


