\name{Z}
\alias{Z}
\docType{data}
\title{
An example of fixed effects (explanatory variables)
}
\description{
An example of fixed effects for 500 samples. The effect consists of three levels, "A", "B", and "C.
}
\details{
Compared with the mean of "A", the mean of "B" is 3 lower and the mean of "C" is 3 higher.
}
\seealso{
  Y, X
}
\examples{
data(sampledata)
dim(Z) #500 samples and 1 effect
unique(Z) #"A", "B", and "C"
}
\keyword{datasets}
