\name{lbsample}
\alias{lbsample}
\title{
Generates gength bias samples from a parametric family.}
\description{ Some famous families are used
}
\usage{
lbsample(n, family, par = list(shape, rate, meanlog, sdlog), censor.vec = rexp(n))
}
\arguments{
  \item{n}{
The sample size.}
  \item{family}{ Some families
}
  \item{par}{Parameters of the family.
}
  \item{censor.vec}{A numeric vector of censoring. When the generated sample is bigger than the values of the vector, the function return 1 for \code{censor}, meaning that sample is right censored.
}
}
\details{Add some details if necessary.
}
\value{
A list containing numeric time of survivals, the censor indicator, and the onset time.}
\references{
Thesis of Pierre-Jerome student.
}
\author{
Pierre-Jerome Bergeron and Vahid Partovi Nia}
\note{
Add some notes about constraints on the parameters.
}

\seealso{
\code{\link{lbfit.nonpar}}, \code{\link{lbfit.par}}.}
\examples{
mydata=lbsample(500,family="exponential",par=list(rate=1))
hist(mydata$data)
}
\keyword{Length bias sampling}
