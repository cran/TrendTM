% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FM_kt.R
\name{FM_kt}
\alias{FM_kt}
\title{It performs the factorization for a fixed rank k and a temporal structure with a fixed tau}
\usage{
FM_kt(
  Data_Series,
  k = 2,
  tau = floor(n/2),
  struct_temp = "none",
  type_soft = "als"
)
}
\arguments{
\item{Data_Series}{the data matrix with d rows and n columns containing the d temporal series with size n.}

\item{k}{the fixed rank of X. Default is 2.}

\item{tau}{the fixed value for tau . Default is \code{floor(n/2)}.}

\item{struct_temp}{a name indicating the temporal structure. Could be \code{none}, \code{periodic} or \code{smooth}. Default is \code{none}.}

\item{type_soft}{the option \code{type} of the function softImpute. Default is \code{als}.}
}
\value{
A list containing
\itemize{
\item \code{M_est} the estimation of M.
\item \code{U_est} the component U of the decomposition of \code{M_est}.
\item \code{V_est} the component V of the decomposition of \code{M_est}.
\item \code{contrast} the Frobenius norm of X-M_est.
}
}
\description{
It performs the factorization for a fixed rank k and a temporal structure with a fixed tau
}
