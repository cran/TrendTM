\name{Data_Series}
\alias{Data_Series}
\docType{data}
\title{Example of data}
\usage{
data(Data_Series)
}
\description{A simulated matrix with a low rank k and with temporal structure based on independent Gaussian.
}
\format{
  A matrix with 30 rows (30 times series) and 100 columns (size of each temporal series).
}
\examples{
library(TrendTM)
data(Data_Series)
head(Data_Series)
TrendTM(Data_Series,k_max=3)
}
\keyword{datasets}
