% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findRowMaxInd.R
\name{findRowMaxInd}
\alias{findRowMaxInd}
\title{Find \eqn{\Sigma} Maximum Indices.}
\usage{
findRowMaxInd(i, max_val, max_ind, row_i, delta, se_est)
}
\arguments{
\item{i}{the row index}

\item{max_val}{the maximal absolute value of row \eqn{i} of the covariance/correlation matrix}

\item{max_ind}{the first index in row \eqn{i} at which \code{max_val} is achieved}

\item{row_i}{a row of \code{abs_sigma}}

\item{delta}{\eqn{\delta}, a numeric constant}

\item{se_est}{vector of standard deviations of features (columns of \code{x})}
}
\value{
a vector of indices
}
\description{
Calculate indices of each row such that the absolute values of these indices
are within \eqn{2\delta} of the maximal absolute value \code{max_val} of this row.
This is an implementation of step 4 of Algorithm 1 in Bing et al. (2020).
}
