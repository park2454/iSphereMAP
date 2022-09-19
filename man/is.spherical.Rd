\name{is.spherical}
\alias{is.spherical}
\title{Check if observations are of unit length}
\description{
  The function checks whether the given matrix has rows with unit length
}
\usage{
  is.spherical(X, tol=1e-6)
}
\arguments{
  \item{X}{matrix}
  \item{tol}{tolerance for numerical errors}
}
\value{
  logical value whether the input matrix is spherical.
}
