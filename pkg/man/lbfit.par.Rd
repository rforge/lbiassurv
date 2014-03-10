\name{lbfit.par}
\alias{lbfit.par}
\title{
Fits length bias corrected survival curve, using parametric method.
}
\description{
The function arguments are time vector, and censoring indicator vector to compute the maximum liklihood estimator of a parametric model.  
}
\usage{
lbfit.par(time, censor, family, initial = list(shape, rate, meanlog, sdlog))
}
\arguments{
  \item{time}{
A vector of positive integers, the survival time.
}
  \item{censor}{
The survival indicator, 0 if censored, 1 for event.
}
  \item{family}{
A string indicating the parametric family to fit. Options are \code{weibull}, \code{gamma}, \code{exponential}, \code{lognormal}
, and \code{loglogistic}.}
  \item{initial}{A list of initial parameters for optimization, \code{meanlog}, and \code{sdlog} for lognormal and \code{shape} and \code{rate} for the rest.
}
}
\details{
This function essentially returns the maximization of likelihood function according to the density function which is indicated by a user in the \code{family} item. The method uses the \code{optim} function with "BFGS" option.
Every element in the \code{time} vector indicates the time that event occurs for one individual, and if the element in the \code{censor} vector is 1 indicates the individual is censored.
}
\value{
\item{shape.fit}{
The estimated value for the shape parameter. 
}
\item{rate.fit}{
The estimated value for the rate parameter.
}
\item{scale.fit}{
If \code{family} is \code{gamma}, \code{scale.fit}=\code{1/rate.fit}.
}
\item{meanlog.fit}{
If \code{family} is \code{lognormal}, it is the estimated mean of the logarithm of a lognormal random variable.
}
\item{sdlog.fit}{
If \code{family} is \code{lognormal}, it is the estimated standard error of the logarithm of a lognormal random variable.
}
\item{shape.sd}{
The estimated standard error for \code{shape} parameter.
}
\item{rate.sd}{
The estimated standard error for \code{rate} parameter.
}
\item{scale.sd}{
The estimated standard error for \code{scale} parameter.
}
\item{meanlog.sd}{
The estimated standard error for \code{meanlog} parameter.
}
\item{sdlog.sd}{
The estimated standard error for \code{sdlog} parameter.
}

\item{loglike.maximum}{
The value of maximization of the likelihood function.
}
\item{hessian}{
A symmetric matrix giving an estimation of the minus Hessian at the solution.
}
\item{code}{
An integer code, 0 indicates successful completion. Other integers indicate possible error codes. See \code{convergence} in \code{\link{optim}} function.
}
  \item{iteration}{
  	A two-element integer vector giving the number of calls related to the optimizing function and its gradient. See \code{\link{optim}} function, counts value return. 
	}
\item{model}{
The model for the fitting algorithm, the same as the  \code{family} parameter.
}
  
  }
	

\references{
Jaime Younger. (2011) Goodness-of-Fit for Length-Biased Survival Data with Right-Censoring.

}
\author{
Pierre-Jerome Bergeron and Vahid Partovi Nia.}

\seealso{
\code{\link{lbfit.nonpar}},\code{\link{lbsample}}.
}
\examples{
mydata=lbsample(500,family="exponential",par=list(rate=1))
lbfit.par(time=mydata$time,censor=mydata$censor,family="exponential",initial=list(shape=3,rate=3))
}
\keyword{length biased sampling}
\keyword{survival curve estimation}