\name{mfpca}
\alias{mfpca}
\title{Multivariate functional pca}
\usage{
  mfpca(fd, nharm, tik = numeric(0))
}
\arguments{
  \item{fd}{in the univariate case fd is an object from a
  class fd. Otherwise in the multivariate case fd is a list
  of fd object (fd=list(fd1,fd2,..)).}

  \item{nharm}{number of harmonics or principal component
  to be retain.}

  \item{tik}{the weights of the functional pca which
  corresponds to the weights of the curves. If don't given,
  then we will run a classic functional pca (without
  weighting the curves).}
}
\value{
  When univarite functional data, the function are
  returning an object of calss "pca.fd", When multivariate
  a list of "pca.fd" object by dimension. The "pca.fd"
  class contains the folowing parameter: harmonics:
  functional data object storing the eigen function values:
  the eigenvalues varprop: the normalized eigenvalues
  (eigenvalues divide by their sum) scores: the scores
  matrix meanfd: the mean of the functional data object
}
\description{
  This function will run a weighted functional pca in the
  two cases of uni, and multivariate cases. If the
  observations (the curves) are given with weights, set up
  the parameter tik.
}
\examples{
data(growth)
data=cbind(matrix(growth$hgtm,31,39),matrix(growth$hgtf,31,54));
t=growth$age;
splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=4);
fd <- Data2fd(data, argvals=t, basisobj=splines);
pca=mfpca(fd,nharm=2)
summary(pca)
}

