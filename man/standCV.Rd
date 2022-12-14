% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standCV.R
\name{standCV}
\alias{standCV}
\title{Standardization for k-fold cross-validation.}
\usage{
standCV(valid_x, train_x, valid_y = NULL, train_y = NULL)
}
\arguments{
\item{valid_x}{the validation set made from the data matrix}

\item{train_x}{the training set made from the data matrix}

\item{valid_y}{the validation set made from the response vector}

\item{train_y}{the training set made from the response vector}
}
\description{
Standardize the training and validation sets for cross-validation.
}
