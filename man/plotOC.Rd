\name{plotOC}
\alias{plotOC}
\title{plot Original Curves}
\usage{
  plotOC(time = 1:nrow(curves), curves, xlab = "time",
    ylab = "value", main = "Original curves")
}
\arguments{
  \item{time}{a vector containing the observation time for
  the curves. If absent the time param will be set at the
  vector 1:nrow(curves)}

  \item{curves}{the observations matrix. Each column of
  this matrix corresponds to one observed curve, and
  contains the value of the curve at discrete time points.}

  \item{xlab}{label of the horizontal axis}

  \item{ylab}{label of the vertical axis}

  \item{main}{the title of the graphic}
}
\description{
  This function plots the observed curves, before any
  action of smoothing or interpolation.
}
\examples{
data(growth)
curves=matrix(data=cbind(growth$hgtm,growth$hgtf),ncol=93)
time=growth$age
plotOC(time,curves)
}

