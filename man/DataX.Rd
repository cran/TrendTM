\name{X}
\alias{X}
\docType{data}
\title{Example of data}
\usage{
data(DataX)
}
\description{A simulated matrix with a low rank k and with temporal structure based on independent Gaussian.
}
\format{
  A matrix with 30 rows (30 times series) and 100 columns (size of each temporal series).
}
\examples{
library(TrendTM)
data(DataX)
head(X)
k.max=3
tau.max=dim(X)[2]
TrendTM(X,k.max=k.max)
}
\keyword{datasets}
