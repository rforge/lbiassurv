\name{lbfit.nonpar}
\alias{lbfit.nonpar}
\title{
Fits a Vardi corrected survival curve, a nonparametric survival cuve estimation.
}
\description{
An extended version of the Kaplan-Miere estimator.
}
\usage{
lbfit.nonpar(time, censor, boot = FALSE, 
boot.control = list(quantile = TRUE, use.median = FALSE, 
confidence.level = 0.95, iter = 1000), 
fit.control = list(iter = 500, tol = 1e-06))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
The observed time as a vector. 
}
  \item{censor}{
The survival indicator, 1 if censored, 0 otherwise.
}
  \item{boot}{
logical, for calculating bootstrap confidence bounds.
}

  \item{boot.control}{
Bootstrap control parameters, as a list.
}

  \item{fit.control}{
Nonparametric fit control parameters.
}
}

\details{
A detailed description, or refer the reader to a work.
}
\value{
The result of the function is an \code{lbsurvfit} object, including step functions for the nonparametric survival curve fit and its bootstrap confidence bounds.}
\references{
Vardi Paper}
\author{
V. P. Nia and P. J. Bergeron}
\note{
Bootstrapping is required if confidence bounds are needed.
}


\seealso{
\code{\link{lbfit.par}}, \code{\link{lbsample}}
}
\examples{
mydata=lbsample(20,family="exponential",par=list(rate=1))
noboot=lbfit.nonpar(time=mydata$time,censor=mydata$censor)
plot(noboot$survfun)
withboot=lbfit.nonpar(time=mydata$time,censor=mydata$censor,boot=TRUE)
x=seq(0,max(mydata$time)+1,length=500)
plot(x,withboot$survfun(x),type="l",col="blue",ylim=c(0,1))
points(x,withboot$lowerfun(x),type="l",col="red")
points(x,withboot$upperfun(x),type="l",col="red")
}

\keyword{Kaplan-Miere estimator}
\keyword{length-bias correction}
