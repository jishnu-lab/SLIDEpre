% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshSigma.R
\name{threshSigma}
\alias{threshSigma}
\title{FDR Thresholding of \eqn{\hat{Sigma}}.}
\usage{
threshSigma(x, sigma, thresh)
}
\arguments{
\item{x}{a data matrix of dimensions \eqn{n \times p}}

\item{sigma}{a sample correlation matrix of dimensions \eqn{p \times p}}

\item{thresh}{a numeric constant}
}
\value{
a matrix of dimensions \eqn{p \times p}
}
\description{
Threshold \eqn{\hat{\Sigma}} according to false discovery rate threshold specification \code{thresh}.
For each entry in \eqn{\hat{\Sigma}}, calculate the associated \eqn{t}-test statistic with degrees of freedom \eqn{n-2}:
\deqn{t = \frac{r\sqrt{n-2}}{\sqrt{1-r^2}}}
Then find the two-tailed p-value for this according to the \eqn{t}-distribution, correct for
multiple testing, and threshold \eqn{\hat{\Sigma}} by setting any entry with an associated p-value
less than or equal to \code{thresh} to 0.
}
