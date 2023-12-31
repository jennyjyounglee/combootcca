% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcs.R
\name{cancor_cov}
\alias{cancor_cov}
\title{Canonical Correlation Analysis for Covariance Matrix}
\usage{
cancor_cov(Sigma, px, align = cca_align_nil, ref)
}
\arguments{
\item{Sigma}{A PSD matrix given the covariance (sample or population)}

\item{px}{The number of variables in the x dataset}

\item{align}{A function to perform post-processing on the estimated
coefficients to render the solution well-identified. By default, this uses
cca_align_posdiag, which ensures that the diagonal of xcoef is
non-negative.}

\item{ref}{Passed through to alignment}
}
\value{
A list containing the following components:

cor: correlations

xcoef: estimated coefficients for the x variables

ycoef: estimated coefficients for the y variables
}
\description{
Compute the Canonical Correlation Analysis solution based on a
covariance matrix
}
\details{
The Canonical Correlation Analysis problem can be solved in a
variety of ways. The approach taken by the cancor function in R
performs QR decomposition on the data matrices and proceeds from
there. However, it is also possible to find the canonical
correlations and associated weights using only the covariance
matrix. This is especially useful if one wishes to solve the
"population" version of CCA for a distribution with population
covariance Sigma. This algorithm results in canonical variates
that have unit empirical variance.
}
\author{
Daniel Kessler
}
