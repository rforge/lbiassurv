\name{vardifit}
\alias{vardifit}
\title{
Implements Vardi correction to the survival function estimation.
}
\description{
This is an EM algorithm type suitable correcting survival curve for the length-biased data.
}
\usage{
vardifit(tt, d, maxiter = 500, tol = 1e-06)
}
\arguments{
  \item{tt}{
what is tt?}
  \item{d}{
what is d?}
  \item{maxiter}{
maximum iteration for the EM algorithm.}
  \item{tol}{
Tolerance for the convergence of the EM algorithm.}
}
\details{
Add details here.
}
\value{
\item{comp1}{Description of the component in the 'list' of the output}
\item{comp2}{Description of the component in the 'list' of the output}
}
\references{
Put the reference here, for instace Vardi paper.
}
\author{
Pierre-Jerome, Vahid.
}
\note{
Add notes here.}
\seealso{
\code{\link{lbiasgenerate}}.}
\examples{
## Generates gamma length-biased data
#survdata <- lbiasgenerate(size=1000)

survdata=lbiasgenerate(100)
bof=vardifit(survdata$data,survdata$censor)
plot(c(0,bof$times),c(1,1-cumsum(bof$ubpvec)),type="s",col="blue")
xx=seq(0,7,by=.01)
lines(xx,dexp(xx),col="red")
xx=seq(0,7,by=.01)
}
\keyword{EM algorithm}
\keyword{Vardi correction}