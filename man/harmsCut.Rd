\name{harmsCut}
\alias{harmsCut}
\title{Separates the matrices of the coefficients of harmonics}
\usage{
  harmsCut(harms, nbasis)
}
\arguments{
  \item{harms}{the coefficients matrix, see the description
  for more details}

  \item{nbasis}{vector containing the number of basis for
  each dimension}
}
\value{
  a list of length dimension of the data which contain the
  harmonic coefficients matrix for each dimension
}
\description{
  This function is used in mfpca to cut the harmonic
  coefficient matrix. In fact mfpca call a c++ mfpca which
  return (among others) the matrix of the coefficients of
  the harmonics in the basis expantion. But in the
  multivariate case the coefficient matrix for all
  dimension are store in the same matrix, so we use this
  function to store each dimension in the specific matrix.
}

