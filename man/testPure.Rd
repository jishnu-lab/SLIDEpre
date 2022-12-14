% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testPure.R
\name{testPure}
\alias{testPure}
\title{Test Node Purity.}
\usage{
testPure(sigma_row, row_ind, s_i, max_vals, max_inds, delta, se_est)
}
\arguments{
\item{sigma_row}{a vector of dimension \eqn{p}}

\item{row_ind}{the index of the row}

\item{s_i}{a vector of indices that are within \eqn{2\delta} of the maximum value of \code{sigma_row}}

\item{max_vals}{a vector of the largest absolute values of each of the rows in \code{sigma}}

\item{max_inds}{a vector of the first index at which each value in Ms is achieved}

\item{delta}{\eqn{\delta}, a numerical constant}

\item{se_est}{a vector of estimates of the standard deviations of the rows of the data matrix, \eqn{x}}
}
\value{
TRUE or FALSE
}
\description{
For given row, check if it is a pure node by iteratively checking the nodes
in \code{Si}. Return TRUE if the given row corresponds to a pure variable.
This is an implementation of steps 6-9 of Algorithm 1 in Bing et al. (2020).
}
