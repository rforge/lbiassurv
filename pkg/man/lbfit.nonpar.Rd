\name{lbfit.nonpar}
\alias{lbfit.nonpar}
\title{
Fits a Vardi corrected survival curve, a nonparametric survival cuve estimation.
}
\description{
An extended version of the Kaplan-Miere estimator.
}
\usage{
lbfit.nonpar(time, censor, initial = list(maxiter = 500, tol = 1e-06))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
The observed time as a vector. 
}
  \item{censor}{
The survival indicator, 1 if censored, 0 otherwise.
}
  \item{initial}{
Some initialization for the nonparametric fit, as a list.}
}
\details{
A detailed description, or refer the reader to a work.
}
\value{
The result of the function is an \code{lbsurvfit} object, including a step function for the nonparametric survival curve fit.}
\references{
Vardi Paper}
\author{
V. P. Nia and P. J. Bergeron}
\note{
Standard errors not provided. Write about confidence bounds.  
}


\seealso{
\code{\link{lbfit.par}}, \code{\link{lbsample}}
}
\examples{
mydata=lbsample(500,family="exponential",par=list(rate=1))
lbfit.nonpar(time=mydata$data,censor=mydata$censor)
}

\keyword{Kaplan-Miere estimator}
\keyword{length-bias correction}
