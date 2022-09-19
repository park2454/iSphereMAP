iSphereMAP = function(X, Y, group_size, lambda=NULL, k_fold=5, rescale=FALSE, verbose=TRUE){
  stopifnot(is.numeric(X))
  stopifnot(!any(is.na(X)))
  stopifnot(dim(X) == dim(Y))
  if(!is.spherical(X) | !is.spherical(Y)){
    if(rescale==TRUE){
      X = as.spherical(X)
      Y = as.spherical(Y)
    } else {
      stop("Input is not spherical")
    }
  }
  n = nrow(X)
  p = ncol(X)
  stopifnot(sum(group_size)==n)
  stopifnot(max(group_size)<= p)
  l2norm = function(x) x^2 |> sum() |> sqrt()
  maxcos = function(x) max(x)/l2norm(x)
  G = length(group_size)
  end = cumsum(group_size)
  start = c(0, end[-G]) + 1
  upper = 1- sqrt(1/2)
  if(is.null(lambda)){
    lambda = seq(-16, log(upper), length=30) |> exp()
  } else if(length(lambda)==1 & is.integer(lambda) & lambda > 0){
    if(lambda==1){
      lambda = upper
    } else {
    lambda = seq(-16, log(upper), length=lambda) |> exp()
    }
  } else {
    lambda = c(lambda[lambda > 0 & lambda < upper], upper )
  }
  if(length(lambda)==0){
    stop("lambda should be an integer or a vector with values between 0 and 1-1/sqrt(2)")
  }
  cv_module = function(lambda){
    shuffle = sample(1:p,p)
    if(p<k_fold){
      K = p
      val_mat = matrix(shuffle,nrow=K)
    } else {
      K = k_fold
      val_mat = matrix(shuffle,nrow=K)
    }
    fold = function(k){
      val_ind = val_mat[k,]
      X_train = X[,-val_ind]
      X_val = X[,val_ind]
      Y_train = Y[,-val_ind]
      Y_val = Y[,val_ind]
      ## Step 1: orthogonal Procustes problem X^TY
      res = svd(t(X_train) %*% Y_train)
      W = res$u %*% t(res$v)
      res2 = svd(t(X_val) %*% Y_val)
      W2 = res2$u %*% t(res2$v)
      ## Step 2: Hard-thresholding
      step2 = function(arg){
        start = arg[1]
        end = arg[2]
        subY = Y_train[start:end,]
        subX = X_train[start:end,]
        fit = lm(t(subY)~t(subX %*% W)+0)
        P = fit$coefficients |> t() |> Matrix()
        colnames(P) = NULL
        tilde_j = apply(P,1,which.max)
        perm = as(tilde_j, "indMatrix")
        if(ncol(perm)<nrow(perm)){
          perm = cbind(perm, Matrix(0, nrow=nrow(perm), ncol = nrow(perm)-ncol(perm)) )
        }
        beta = 1 - apply(P,1,max)/apply(P,1,l2norm)
        ind = beta > lambda
        P[ind,] = P[ind,] / apply(P[ind,]%*%subX,1,l2norm)
        P[!ind,] = perm[!ind,]
        (Y_val[start:end,] - P %*% X_val[start:end,] %*% W2)^2 |> sum()
      }
      cbind(start,end) |> apply(1,step2) |> sum()
    }
    lapply(1:K,fold) |> unlist() |> mean()
  }
  if(length(lambda)>1){
    if(verbose==TRUE){
      cv_error = pblapply(lambda,cv_module) |> unlist()
    } else {
      cv_error = lapply(lambda,cv_module) |> unlist()
    }
    lambda_cv = lambda[which.min(cv_error)]
  } else {
    cv_error = NULL
    lambda_cv = lambda
  }
  ## Step 1: orthogonal Procustes problem X^TY
  res = svd(t(X) %*% Y)
  W = res$u %*% t(res$v)
  ## Step 2: Hard-thresholding
  step2 = function(arg){
    start = arg[1]
    end = arg[2]
    subY = Y[start:end,]
    subX = X[start:end,]
    fit = lm(t(subY)~t(subX %*% W)+0)
    P = fit$coefficients |> t() |> Matrix()
    colnames(P) = NULL
    tilde_j = apply(P,1,which.max)
    perm = as(tilde_j, "indMatrix")
    if(ncol(perm)<nrow(perm)){
      perm = cbind(perm, Matrix(0, nrow=nrow(perm), ncol = nrow(perm)-ncol(perm)) )
    }
    beta = 1 - apply(P,1,max)/apply(P,1,l2norm)
    ind = beta > lambda_cv
    P[ind,] = P[ind,] / apply(P[ind,]%*%subX,1,l2norm)
    P[!ind,] = perm[!ind,]
    correct = !ind & (tilde_j == seq(length(tilde_j)))
    return(list(P=P, correct = correct))
  }
  P_list = cbind(start,end) |> apply(1,step2)
  P_large = Diagonal(n)
  pb <- txtProgressBar(min = 0,
                       max = G,
                       style = 3,
                       width = 50,
                       char = "=")
  for(g in 1:G){
    P_large[start[g]:end[g],start[g]:end[g]] = P_list[[g]][[1]]
    if(verbose==TRUE) setTxtProgressBar(pb, g)
  }
  close(pb)
  correct = lapply(1:G, function(g) P_list[[g]][[2]]) |> unlist()
  ## Step 3: orthogonal Procustes problem X^TY
  res = svd(t(X[correct,]) %*% Y[correct,])
  W = res$u %*% t(res$v)
  output = list(w = W, p = P_large, lambda = lambda_cv,
                correct=correct, lambda_grid = lambda, cv_error = cv_error)
  return(output)
}
