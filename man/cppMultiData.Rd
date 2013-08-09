\name{cppMultiData}
\alias{cppMultiData}
\title{The C++ code, of this package take to run from the data two main information.
the coefficient in the basis expansion of the functional data, and the inner product between theses basis.
cppMultidata made this task in the multivariate case.}
\usage{
  cppMultiData(mfd)
}
\arguments{
  \item{mfd}{a list containing all dimension of the data}
}
\value{
  fdData a list that containing the concatenation of the
  transpose of the coefficients matrix and a block diagonal
  matrix where each block is the inner product between the
  basis functions.
}
\description{
  The C++ code, of this package take to run from the data
  two main information. the coefficient in the basis
  expansion of the functional data, and the inner product
  between theses basis. cppMultidata made this task in the
  multivariate case.
}

