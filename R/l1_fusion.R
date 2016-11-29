#require("irlba")
#require("rARPACK")

#' Calculate maximal eigenvalue for big matrices using singular
#' value decomposition
#'
#' @param X Matrix to be evaluated (can be a bigmatrix)
#' @param method One of 'irlba' or 'rARPACK'
#'
#' @return The maximal eigenvalue 
#' @export
#'
#' @examples
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

#' Fused lasso optimisation with the proximal-gradient method (Chen et al. 2010)
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
#' @export
#'
#' @examples
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
    
    # cat(L_U, '\n')
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
        
        if (improvement < tol) 
            break
        
        B.old = B.new
        
        
    }
    
    if (i >= num.it) 
        warning("Reached max iterations without convergence.")
    
    return(B.new)
}


# Fused lasso optimisation with proximal-gradient method -- fixed constraints
fusedLassoProximal <- function(X, Y, groups, lambda, gamma, G, mu = 1e-04, tol = 1e-06, num.it = 1000, lam.max = NULL, c.flag = FALSE, intercept = TRUE, penalty.factors = rep(1, 
    dim(X)[2])) {
    
    group.names = sort(unique(groups))
    samp.sizes = table(groups)[group.names]
    
    k = length(group.names)
    p = dim(X)[2]
    
    if (intercept) 
        X = cbind(X, 1)
    
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
        
        X = NULL
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
        c.flag, intercept = intercept, penalty.factors = penalty.factors, X.list = X.list))
}
