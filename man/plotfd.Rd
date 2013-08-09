\name{plotfd}
\alias{plotfd}
\title{plot a functional data object}
\usage{
  plotfd(fd, col = c(1:nrow(fd$coefs)), xlab = "time",
    ylab = "value", main = "Functional data curves")
}
\arguments{
  \item{fd}{a functional data object}

  \item{col}{the color vector.}

  \item{xlab}{label of the horizontal axis}

  \item{ylab}{label of the vertical axis}

  \item{main}{the title of the graphic}
}
\description{
  This function plots a functional data object (after
  smoothing or interpolation). If you want to color the
  curves according to a cluster's membership, please
  specify the parmetere col. Note: this function works only
  for univariate functional data.
}
\examples{
data(growth)
data=cbind(matrix(growth$hgtm,31,39),matrix(growth$hgtf,31,54));
t=growth$age;
splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=4);
fd <- Data2fd(data, argvals=t, basisobj=splines);
cls=c(rep(1,39),rep(2,54)) #there is 39 boys and 54 girls in this data
plotfd(fd,col=cls)
}

