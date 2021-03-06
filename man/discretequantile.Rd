% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdp_basic_functions.R
\name{discretequantile}
\alias{discretequantile}
\title{quantile function of a discrete random variable Y that takes on values y with probabilities py}
\usage{
discretequantile(p, y, py, ordered = FALSE)
}
\arguments{
\item{p}{A number in (0,1) representing the desired quantile (p=0.5 correspons to the median).}

\item{y}{A vector representing the rank of a discrete random variable.}

\item{py}{A vector of with elements in (0,1) that add up to 1, representing the probabilities that Y=y}

\item{ordered}{Logical; indicates whether or not the vector y is sorted in ascending order}

\item{A}{list with elements x (the matrix of covariates), y (outcomes) and a (the missingness indicator)}

\item{Logical;}{if TRUE, the vector y is assumed ordered, if FALSE it is not.}
}
\value{
The \code{p}-quantile of Y
}
\description{
quantile function of a discrete random variable Y that takes on values y with probabilities py
}
