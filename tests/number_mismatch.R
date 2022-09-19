library(iSphereMAP)

p = 300
alpha = 0.8
n = 2000
K= 100

generator = function(n,p,alpha,K){
  X = rnorm(n*p) |> matrix(n,p) |> as.spherical()
  W = rnorm(p^2) |> matrix(p,p)
  W = svd(W)$u
  nmis = sample(1:(n/K),K, replace=T)
  nmis = round(nmis/sum(nmis)*(n^alpha))
  selector = function(k){
    ind = 1:(n/K)
    mis = sample(ind,nmis[k])
    ind[mis] = sample(mis,nmis[k])
    ind = ind + (k-1)*(n/K)
    return(ind)
  }
  ind = sapply(1:K,selector) |> as.vector()
  Z = X[ind,]
  Y = Z%*%W
  return(list(x=X, y=Y, z=Z, nmis=sum(ind != 1:n), ind=ind))
}
gen = generator(n,p,alpha,K)
fit = iSphereMAP(gen$x, gen$y, rep(n/K,K))
fit$nmis == gen$nmis
