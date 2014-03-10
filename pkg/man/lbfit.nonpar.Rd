\name{lbfit.nonpar}
\alias{lbfit.nonpar}
\title{
Fits length bias corrected survival curve, using nonparametric method estimation.
}
\description{
This function implements length-bias correction to the nonparametricestimator.
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
The survival indicator, 0 if censored, 1 for event.
}
  \item{boot}{
logical, for calculating bootstrap confidence bounds. By default is \code{FALSE}.
}

  \item{boot.control}{
Bootstrap control parameters, as a list. See 'Details'.
}

  \item{fit.control}{
Nonparametric fit control parameters, as a list. See 'Details'.
}
}

\details{
If \code{boot=TRUE}, the survival function is reported with two additional functions, one for lower bound, and another for the upper bound. While \code{boot=FALSE}, no standard error and no confidence bounds are produced. The \code{boot.control} items and \code{fit.control} items can be changed, if is necessary. The \code{quantile} gives the quantile-based confidence bounds, and by default is \code{TRUE}. If \code{quantile=FALSE} a normal confidence bound is produced. The \code{use.median} replaces the mean with median as the center of the confidence bounds, when \code{quantile=FALSE}. The \code{confidence.level} is the level of confidence for one bootstrap interval. The \code{iter} is the number of iterations for computing the bootstrap confidence bounds. The \code{tol} is the convergence tolerance.
}
\value{
\item{par}{A matrix, indicates survival estimate for each observed  time parameter. 
 }
\item{survfun}{
The survival function as a step function. 
}
\item{iterations}{
  An integer value giving the number of calls.
}
\item{conv}{
A logical value, it is \code{TRUE} if the function converges, and if not, the value is \code{FALSE}.
}
\item{method}{Equals \code{nonpar}. This value will be used in constructing a proper survival function plot.
}
\item{lowerfun}{
A step function for the lower bound of confidence interval for the nonparametric estimator.
}
\item{upperfun}{A step function for the upper bound of confidence interval for the nonparametric estimator.
}
\item{attr}{
The class of the object.
}
}
\references{
Vardi, Y. (1989). Multiplicative censoring, renewal processes, deconvolution and decreasing density: Nonparametric estimation.  \emph{Biometrika}, \bold{76} (4), 751--761.}
\author{
Pierre-Jerome Bergeron and Vahid Partovi Nia.}
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