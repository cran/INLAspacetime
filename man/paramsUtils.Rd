% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paramsUtils.R
\name{paramsUtils}
\alias{paramsUtils}
\alias{gsConstant}
\alias{gammas2params}
\alias{params2gammas}
\title{Funtions to help converting from/to user/internal parametrization.}
\usage{
gsConstant(lgammas, alpha, smanifold)

gammas2params(lgammas, alpha.t, alpha.s, alpha.e, smanifold = "R2")

params2gammas(lparams, alpha.t, alpha.s, alpha.e, smanifold = "R2")
}
\arguments{
\item{lgammas}{the SPDE parameters log(gamma.s, gamma.t, gamma.e)}

\item{alpha}{the resulting spatial order.
Values could be "S1", "S2", "R1", "R2" and "R3".}

\item{smanifold}{spatial domain manifold.}

\item{alpha.t}{temporal order of the SPDE}

\item{alpha.s}{spatial order of the spatial differential operator
in the non-separable part.}

\item{alpha.e}{spatial order of the spatial differential operator
in the separable part.}

\item{lparams}{log(spatial range, temporal range, sigma)}
}
\value{
the part of sigma due to spatial constant and gamma.s

log(spatial range, temporal range, sigma)

log(gamma.s, gamma.t, gamma.e)
}
\description{
Funtions to help converting from/to user/internal parametrization.
}
\examples{
gammas2params(log(c(0, 0, 0)), 1, 2, 1, "R2")
params2gammas(log(c(1, 1, 1)), 1, 2, 1, "R2")
}
