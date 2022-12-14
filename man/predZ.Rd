% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predZ.R
\name{predZ}
\alias{predZ}
\title{Estimate \eqn{Z}.}
\usage{
predZ(x, er_res)
}
\arguments{
\item{x}{data matrix, standardized}

\item{er_res}{the results of a run of \code{plainER()} or \code{priorER()}}
}
\value{
estimates for \eqn{Z}
}
\description{
Additional prediction function using the support of \eqn{\beta}.
}
