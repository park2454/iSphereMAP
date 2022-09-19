\name{iSphereMAP}
\alias{iSphereMAP}
\title{Three-step estimator from the iSphereMAP procedure}
\description{
Three-step estimator from the iSphereMAP procedure. Returns the orthogonal regression map and the permutation matrix simultaneously.
}
\usage{
  iSphereMAP(X, Y, group_size, lambda=NULL, k_fold=5, rescale=FALSE, verbose=TRUE)
}
\arguments{
  \item{X}{predictor matrix, required to be spherical unless \code{rescale=TRUE}}
  \item{Y}{response matrix, required to be spherical unless \code{rescale=TRUE}}
  \item{group_size}{a numeric vector representing the group information of input. The components of the vector represents the number of observations in the corresponding group. Note that the individual group size has to be smaller than the number of columns in X and Y}
  \item{lambda}{a positive integer or a numeric vector. Given a positive integer a grid of values with length equal to that integer will be generated. Users can manually provide a sequence of values as well. Note that the tuning parameter should be between 0 and 1-1/sqrt(2). By default, a grid of length 20 will be generated.}
  \item{k_fold}{positive integer, the number of folds to be used in K-fold cross-validation}
  \item{rescale}{logical; If TRUE, rescales each observations in X and Y to have unit value. If FALSE and either X or Y is not spherical, the function will return an error. FALSE by default.}
  \item{verbose}{logical; If TRUE, a progress bar for cross-validation and mismatch collecting process is shown}
}
\value{
  \item{w}{estimated regression coefficient. an orthogonal matrix}
  \item{p}{estimated permutation matrix. If mismatch is sparse, the permutation matrix will be represented as a sparse matrix}
  \item{lambda}{tne value of the tuning parameter in the estimated model (chosen so that the model has at most \code{max.mismatch}) mismatches}
  \item{mismatch}{estimated number of observations that are mismatch}
}
\source{
  The iSphereMAP function relies on \code{lm}, \code{as} and \code{Matrix} library.
}
\references{
  Shi, X., Li, X., & Cai, T. (2021). Spherical regression under mismatch corruption with application to automated knowledge translation. Journal of the American Statistical Association, 116(536), 1953-1964.
}
\examples{
n=1000
p=50
X = rnorm(n*p) |> matrix(n,p) |> as.spherical()
W = rnorm(p^2) |> matrix(p,p)
W = svd(W)$u
shuffle = 1:n
mismatch = sample(shuffle,n^0.5)
shuffle[mismatch] = sample(mismatch,length(mismatch))
Z = X[shuffle,]
Y = Z%*%W
fit = iSphereMAP(X,Y,group_size = rep(20,n/20))
fit$lambda
A = sum(shuffle != 1:n) # number of mismatch in the data
B = fit$mismatch # estimated number of mismatch in the data
A == B; A
}
