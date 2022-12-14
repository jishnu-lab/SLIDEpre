% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getOptPerm.R
\name{getOptPerm}
\alias{getOptPerm}
\title{Find Best Column Permutation Matrices.}
\usage{
getOptPerm(A, B)
}
\arguments{
\item{A}{a matrix}

\item{B}{a matrix of same dimensions as \eqn{A}}
}
\value{
a list of column and sign permutation matrices
}
\description{
For the given pure node set, find the best column permutation matrix
\eqn{A_{perm}} of \eqn{A} such that \eqn{|A_{perm} - B|_F} is minimal for the target \eqn{B}.
}
