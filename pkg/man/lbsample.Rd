\name{lbsample}
\alias{lbsample}
\title{
 Generates length bias samples from a parametric family.}
\description{ Generates length bias samples from a parametric family. 
}
\usage{
lbsample(n, family, par = list(shape, rate, meanlog, sdlog), censor.vec = rexp(n))
}
\arguments{
  \item{n}{
An integer indicates the number of the desired samples.}
  \item{family}{ A string. Options are \code{weibull}, \code{gamma}, \code{exponential}, \code{lognormal}
, and \code{loglogistic}.}
  
\item{par}{A list of parameters related to the chosen \code{family},  \code{meanlog}, and \code{sdlog} for \code{lognormal}, and \code{shape} and \code{rate} for the rest.
}
  \item{censor.vec}{An integer vector of the censoring indicator. When the generated sample is bigger than the values of the vector, the function return 1 in the output, meaning that sample is right censored. 
}
}
\details{
This function uses the sampling scheme described in Younger(2011).\cr


{\code{weibull}:   f(x)= k*b*(x*b)^(k-1)*exp(1(x*b)^k)}\cr

x>=0, rate=b>0, shape=k>0.\cr

{\code{exponential}:   f(x)=lambda*exp(-lambda*x)}\cr

x>=0, rate=lambda>0.\cr

{\code{lognormal}:}    f(x)=[1/(x*sigma*sqrt(2*pi))]*exp(-[(ln(x)-mu)^2]/2*sigma^2)\cr

x>0, sdlog=sigma>0, meanlog=mu.\cr

{\code{loglogistic}:    f(x)=[(b/a)*(x/a)^(b-1)]/[1+(x/a)^b]^2\cr

x>=0, rate=1/a>0, shape=b>0.}\cr

  
}


\value{
\item{time}{A vector of \code{n} random length bias samples from the  density function specified in \code{family}, giving the residual lifetime.}
\item{censor}{A vector of \code{n} boolean integers to indicates if the sample is censored or not, 1 for censored samples and 0 for uncensored.}
\item{onset}{
Starting time of the follow up for each individual, sampled from uniform distribution.}
}
\references{
Jaime Younger. (2011) Goodness-of-Fit for Length-Biased Survival Data with Right-Censoring.
}
\author{
Pierre-Jerome Bergeron and Vahid Partovi Nia}
\note{
The length of the \code{censor.vec} must match the number of the desired samples.
}

\seealso{
\code{\link{lbfit.nonpar}}, \code{\link{lbfit.par}}.}
\examples{
mydata=lbsample(500,family="exponential",par=list(rate=1))
hist(mydata$time)
}
\keyword{Length bias sampling}