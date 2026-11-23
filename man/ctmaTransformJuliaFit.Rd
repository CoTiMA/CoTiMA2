\name{ctmaTransformJuliaFit}
\alias{ctmaTransformJuliaFit}
\title{Transform Julia Point Estimates Without Calculating Uncertainty}
\usage{
ctmaTransformJuliaFit(fit, setup, cores = 1L, digits = 4L)
}
\arguments{
\item{fit}{A completed \code{ctmaJuliaFit} object.}
\item{setup}{A matching fitted, \code{fit=FALSE}, or prepared ctsem/CoTiMA
setup object supplying ctsem's transformation model.}
\item{cores}{Number of cores passed to ctsem's transformation code.}
\item{digits}{Number of decimal places used internally by the ctsem summary.}
}
\value{
The updated \code{ctmaJuliaFit} object. Point transformations are stored in
\code{coefficients$transformed}, \code{coefficients$population},
\code{coefficients$tipred}, \code{coefficients$matrices}, and, where
available, \code{coefficients$study}. No covariance matrix, transformed draws,
or uncertainty summary is added.
}
\description{
Applies ctsem's parameter constraints and generated-quantity transformations
to the raw maximum-likelihood estimates stored in a completed
\code{ctmaJuliaFit}. Unlike \code{ctAddJuliaSamples()}, this function does not
require a Hessian or covariance matrix, draw samples, or repeat optimization.
}
\details{
The resulting transformed components are point estimates only. They must not
be interpreted as standard errors, confidence intervals, or evidence that the
optimizer converged.
}
