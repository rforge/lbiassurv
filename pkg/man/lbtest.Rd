\name{lbtest}
\alias{lbtest}
\title{This function implements asymptotic and bootstrap test of length bias.}
\description{This function implements a visual and a formal test. The formal test produces the asymptotic and the bootstrap p-values.}
\usage{
lbtest(time, censor, entry, offset = 0, bootstrap = TRUE, plot = TRUE, alternative = "two.sided", boot.control = list(confidence.level = 0.95, iter = 1000), plot.control = list())
}
\arguments{
  \item{time}{The observed time as a vector.}
  \item{censor}{The survival indicator, 1 if censored, 0 otherwise.}
  \item{entry}{Time to entry.}
  \item{offset}{Offset value.}
  \item{bootstrap}{Logical, generate bootstrap p-values.}
  \item{plot}{Logical, generate a visual test.}
  \item{alternative}{The alternative hypothesis for computing the p-value, two.sided, less, or greater.}
  \item{boot.control}{A list of control parameters for bootstrap.}
  \item{plot.control}{A list of control parameters for the visual test.}
}
\details{The bootstrap is time consuming for large data.}
\value{A list of test statistic, its standard error and the p-value.
}
\references{Addona & Wolfson (2006), Wei (1980), Asgharian, Wolfson & Zhang (2006).}
\author{Pierre-Jerome}
\examples{
mydata=lbsample(1000,family="exponential",par=list(rate=1))

lbtest(time=mydata$time, censor=mydata$censor, entry=mydata$entry, 
offset=0, bootstrap=FALSE, plot=FALSE,alternative="two.sided")

}
\keyword{length-biased stationarity test}
