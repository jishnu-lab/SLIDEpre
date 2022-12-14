% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustSign.R
\name{adjustSign}
\alias{adjustSign}
\title{Adjust Matrix Signs.}
\usage{
adjustSign(sigma, AI)
}
\arguments{
\item{sigma}{a sample correlation matrix of dimensions \eqn{p \times p}}

\item{AI}{a matrix of dimensions \eqn{p \times K}}
}
\value{
a matrix of dimensions \eqn{p \times p}
}
\description{
Perform a sign operation on matrix \eqn{\Sigma} according to the sign of \eqn{A_I}.
Rows corresponding to nodes not in \eqn{A_I} are set to 0.
}
