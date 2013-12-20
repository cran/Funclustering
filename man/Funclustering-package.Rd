\name{Funclustering-package}
\alias{Funclustering-package}
\alias{Funclustering}
\docType{package}
\title{
A package for functional data clustering.
}
\description{
This packages proposes a model-based clustering algorithm for multivariate functional data.
The parametric mixture model, based on the assumption of normality of the principal components resulting from a multivariate functional PCA, is estimated by an EM-like algorithm.
The main advantage of the proposed algorithm is its ability to take into account the dependence among curves. 


}
\details{

This package include a function named funclust which allow to perform functional data clustering. 
This package allows also to perform a functional PCA (function mfpca) for univariate or multivariate functional data, with possibility to define different weights for each observation.  
The output of funclust and mfpca are a list, so one can use summary() to display the results.
We can also plot the original curves with function plotOC. The function plotfd plots the curves after interpolation or smoothing.
See the help of theses functions for more details.  
}

\references{
J.Jacques and C.Preda (2013), Funclust: a curves clustering method using functional random variable density approximation, Neurocomputing, 112, 164-171.

J.Jacques and C.Preda (2013), Model-based clustering of multivariate functional data, Computational Statistics and Data Analysis, in press DOI 10.1016/j.csda.2012.12.004.
}

\examples{
 # Multivariate
# ---------  CanadianWeather (data from the fda package) --------
 CWtime<- 1:365
 CWrange<-c(1,365)
 CWbasis <- create.fourier.basis(CWrange, nbasis=65)
 harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), rangeval=CWrange)
 
 # -- Build the curves --
 temperature=CanadianWeather$dailyAv[,,"Temperature.C"]
 CWfd1 <- smooth.basisPar(
 CWtime, CanadianWeather$dailyAv[,,"Temperature.C"],CWbasis,
 Lfdobj=harmaccelLfd, lambda=1e-2)$fd
 precipitation=CanadianWeather$dailyAv[,,"Precipitation.mm"]
 CWfd2 <- smooth.basisPar(
 CWtime, CanadianWeather$dailyAv[,,"Precipitation.mm"],CWbasis,
 Lfdobj=harmaccelLfd, lambda=1e-2)$fd
 
 # -- the multivariate functional data object --  
 CWfd=list(CWfd1,CWfd2) 
 
 # -- clustering in two class -- 
 res=funclust(CWfd,K=2)
 summary(res)

}
