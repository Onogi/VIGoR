\name{Y}
\alias{Y}
\docType{data}
\title{
An example of response variables.
}
\description{
A vector consisting of 500 samples.
}
\details{
Y mimicked phenotypic values in quantitative genetics. Y was created from X (SNP genotypes) and Z (fixed effect). 10 SNPs in X have additive effects, and 5 out of 10 SNPs have dominance effects further. Y also include environmental noises. The variances of additive, dominance, and noise are approximately 1 : 1 : 1.
}
\seealso{
  X, Z
}
\examples{
data(sampledata)
length(Y) #500 samples
any(is.na(Y)) #FALSE. Y has no missing (but it's allowed).
}
\keyword{datasets}
