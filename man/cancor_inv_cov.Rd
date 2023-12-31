% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcs.R
\name{cancor_inv_cov}
\alias{cancor_inv_cov}
\title{Inverse Map from CCA to Covariance}
\usage{
cancor_inv_cov(xcoef, ycoef, cor)
}
\arguments{
\item{xcoef}{Canonical directions for the x variables (matrix)}

\item{ycoef}{Canonical directions for the y variables (matrix)}

\item{rho}{Canonical correlations (vector)}
}
\value{
A square covariance matrix of size p + q
}
\description{
Given canonical directions and correlations, find an associated covariance matrix
}
\details{
In the same way that cancor_cov maps from a covariance matrix to a
correspond CCA solution, it is possible to map from a CCA solution to a
covariance matrix. However, this map is generally not unique.

Given canonical correlations along with canonical directions for x and y,
respectively, this function will return a covariance matrix Sigma such that
calling cancor_cov on Sigma will return the same parameterization.

Note that when p is not equal to q, there is ambiguity in the specification
of the covariance matrix for the larger random vector. The current approach
involves pseudoinverses, but these present difficulties for downstream
matrix operations. As a fix, we "fortify" the otherwise rank-deficient
matrices by increasing their trailing eigenvalues.
}
\author{
Dan Kessler
}
