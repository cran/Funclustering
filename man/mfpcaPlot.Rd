\name{mfpcaPlot}
\alias{mfpcaPlot}
\title{Plot multivariate functional pca}
\usage{
  mfpcaPlot(pca, grid = c())
}
\arguments{
  \item{pca}{is the result of mfpca. In the univariate case
  mfpcaPlot use the package fda and will be similar to it's
  function "plot.pca.fd". In multivariate functional pca,
  we will make a graphic window for each dimension.}

  \item{grid}{specify how to divide the graphics window.
  grid=c(n,m) divided the widow in to n lines and m
  columns. If user don't specify grid then he must enter
  <Enter> to pass to the next graphic.}
}
\description{
  This function plots the functional pca.
}
\examples{
# Multivariate
# ---------  CanadianWeather (data from the package fda) --------
CWtime<- 1:365
CWrange<-c(1,365)
CWbasis <- create.fourier.basis(CWrange, nbasis=65)
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), rangeval=CWrange)

# -- Build the curves ---
temp=CanadianWeather$dailyAv[,,"Temperature.C"]
CWfd1 <- smooth.basisPar(
CWtime, CanadianWeather$dailyAv[,,"Temperature.C"],CWbasis,
Lfdobj=harmaccelLfd, lambda=1e-2)$fd
precip=CanadianWeather$dailyAv[,,"Precipitation.mm"]
CWfd2 <- smooth.basisPar(
CWtime, CanadianWeather$dailyAv[,,"Precipitation.mm"],CWbasis,
Lfdobj=harmaccelLfd, lambda=1e-2)$fd

CWfd=list(CWfd1,CWfd2)

pca=mfpca(CWfd,nharm=4)
mfpcaPlot(pca,grid=c(2,2))
}

