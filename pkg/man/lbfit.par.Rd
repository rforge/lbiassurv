\name{lbfit.par}
\alias{lbfit.par}
\title{
Fits parametric models
}
\description{
Will be combined with nonpar version soon.
}
\usage{
lbfit.par(time, censor, family, initial = list(shape, rate, meanlog, sdlog))
}
\arguments{
  \item{time}{
The observed time as a vector. 
}
  \item{censor}{
The survival indicator, 1 if censored, 0 otherwise.
}
  \item{family}{
A string indicating the parametric family to fit. Options are \code{weibull}, \code{gamma}, \code{exponential}, \code{lognormal}
, and \code{loglogistic}.}
  \item{initial}{A list of initial parameters for optimization, \code{meanlog}, and \code{sdlog} for lognormal and the \code{shape} and \code{rate} for the rest.
}
}
\details{
See the master thesis of Pjerome student.
}
\value{
The list returned by the \code{optim} function.}
\references{
The thesis or an article here.
}
\author{
V.P. Nia and and P. J. Bergeron.}
\note{
Write why and when you get NaNs.  
}

\seealso{
\code{\link{lbfit.nonpar}},\code{\link{lbsample}}.
}
\examples{
mydata=lbsample(500,family="exponential",par=list(rate=1))
lbfit.par(time=mydata$data,censor=mydata$censor,family="exponential",initial=list(shape=3,rate=3))
}
\keyword{length biased sampling}
\keyword{survival curve estimation}