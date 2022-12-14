% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calFittedSigma.R
\name{calFittedSigma}
\alias{calFittedSigma}
\title{Calculate Fitted Sigma.}
\usage{
calFittedSigma(sigma, delta, max_vals, max_inds, se_est)
}
\arguments{
\item{sigma}{a correlation matrix of dimensions \eqn{p \times p}}

\item{delta}{a threshold parameter}

\item{max_vals}{the calculated maximal values of \eqn{\Sigma} by row}

\item{max_inds}{vector of column indices at which the values in \code{max_vals} are achieved in \code{| sigma |}}

\item{se_est}{the estimated standard errors}
}
\value{
a list containing a vector of the indices of the estimated pure variables and
the fitted value of \eqn{A_I \cdot C \cdot A_I^\top}. returns -1 if only one pure node is identified
}
\description{
Calculate the fitted value of \eqn{A_I \cdot C \cdot A_I^\top} for given
\eqn{\Sigma} and \eqn{\delta}.
}
