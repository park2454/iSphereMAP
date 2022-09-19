is.spherical = function(X, tol=1e-6){
  # check if rows of design matrix is of unit length
  if(!is.numeric(X)){
    message("Input is not numeric")
    return(FALSE)
  }
  if(any(is.na(X))){
    message("Input contains NA")
    return(FALSE)
  }
  l2norm = function(x) x^2 |> sum() |> sqrt()
  if( any(abs(apply(X,1,l2norm)-1)>1e-6) ){
    message("Some observations are not of unit length, use as.spherical()")
    return(FALSE)
  }
  return(TRUE)
}
