as.spherical = function(X){
  # scale observations in the design matrix to unit length
  stopifnot(is.numeric(X))
  stopifnot(!any(is.na(X)))
  l2norm = function(x) x^2 |> sum() |> sqrt()
  X / apply(X,1,l2norm)
}
