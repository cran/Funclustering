\name{funclust}
\alias{funclust}
\title{funclust, clustering multivariate functional data}
\usage{
  funclust(fd, K, thd = 0.05, increaseDimension = FALSE,
    hard = FALSE, fixedDimension = integer(0), nbInit = 5,
    nbIterInit = 20, nbIteration = 200, epsilon = 1e-05)
}
\arguments{
  \item{fd}{in the univariate case fd is an object from a
  class fd of fda package. Otherwise, in the multivariate
  case, fd is a list of fd objects (fd=list(fd1,fd2,..)).}

  \item{K}{the number of clusters.}

  \item{thd}{the threshold in the Cattell scree test used
  to select the numbers of principal components retained
  for the approximation of the probability density of the
  functional data.}

  \item{increaseDimension}{A logical parameter. If FALSE,
  the numbers of principal components are selected at each
  step of the algorithm according to the Cattell scree
  test. If TRUE, only an increase of the numbers of
  principal components are allowed.}

  \item{hard}{A logical parameter. If TRUE, the algorithm
  is randomly initialized with "hard" cluster membership
  (each curves is randomly assigned to one of the
  clusters). if FALSE, "soft" cluster membership are used
  (the cluster membership probabilities are randomly
  choosen in the K-simplex).}

  \item{fixedDimension}{A vector of size K which contains
  the numbers of principal components in the case of fixed
  numbers of principal components.}

  \item{epsilon}{The stoping criterion for the EM-like
  algorithm: the algorithm is stopped if the difference
  between two successive loglikelihood is less than
  epsilon.}

  \item{nbInit}{The number of small-EM used to determine
  the intialization of the main EM-like algorithm.}

  \item{nbIterInit}{The maximum number of iterations for
  each small-EM.}

  \item{nbIteration}{The maximum number of iteration in the
  main EM-like algorithm.}
}
\value{
  a list containing: tik: conditional probabilities of
  cluster membership for each observation cls: the
  estimated partition, proportions: the mixing proportions,
  loglikelihood: the final value of the approximated
  likelihood, loglikTotal: the values of the approximated
  likelihood along the EM-like algorithm (not necessarily
  increasing), aic, bic, icl: model selection criteria,
  dimensions: number of principal components, used in the
  functional data probability density approximation,
  selected for each clusters, dimTotal: values of the
  number of principal components along the EM-like
  algorithm, V: principal components variances per cluster
}
\description{
  This function run clustering algorithm for multivariate
  functional data.
}
\details{
  There is multiple kind of running the function funclust.
  The first one is to run the function with fixed
  dimensions among all iteration of the algorithm.
  (parameter fixedDimenstion). fixedDimension must be an
  integer vector of size K (number of cluster). If the user
  gives a not integer value in fixedDimension, then this
  value will be convert automatically to an integer one for
  examplethe, the algorithm will run with a dimension 2
  instead of 2.2 given by user. The second one is to run it
  with a dimensions variying according to the results of
  the scree test (parameter increaseDimension). If
  increaseDimension = true, then the dimensions will be
  constraint to only increase between to consecutuves
  iterations of the algorithm. else the values of the
  dimensions will be the results of the scree test.
}
\examples{
data(growth)
data=cbind(matrix(growth$hgtm,31,39),matrix(growth$hgtf,31,54));
t=growth$age;
splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=4);
fd <- Data2fd(data, argvals=t, basisobj=splines);
# with varying dimensions (according to the results of the scree test)
res=funclust(fd,K=2)
summary(res)

# with fixed dimensions
res=funclust(fd,K=2,fixedDimension=c(2,2))


# Multivariate (deactivated by default to comply with CRAN policies)
# ---------  CanadianWeather (data from the package fda) --------
# CWtime<- 1:365
# CWrange<-c(1,365)
# CWbasis <- create.fourier.basis(CWrange, nbasis=65)
# harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), rangeval=CWrange)

# -- Build the curves ---
# temp=CanadianWeather$dailyAv[,,"Temperature.C"]
# CWfd1 <- smooth.basisPar(
# CWtime, CanadianWeather$dailyAv[,,"Temperature.C"],CWbasis,
# Lfdobj=harmaccelLfd, lambda=1e-2)$fd
# precip=CanadianWeather$dailyAv[,,"Precipitation.mm"]
# CWfd2 <- smooth.basisPar(
# CWtime, CanadianWeather$dailyAv[,,"Precipitation.mm"],CWbasis,
# Lfdobj=harmaccelLfd, lambda=1e-2)$fd

# CWfd=list(CWfd1,CWfd2)

# res=funclust(CWfd,K=2)
}

