\name{ctmaStudyEstimates}
\alias{ctmaStudyEstimates}
\title{Reconstruct native Julia study DRIFT estimates and joint covariance}
\usage{ctmaStudyEstimates(fit, cores = 1L, derivativeStep = 1e-5)}
\arguments{
\item{fit}{A native ctsem Julia ML fit or pooled CoTiMA2 fit with saturated study contrasts and no hierarchical or person-level random effects.}
\item{cores}{Positive integer number of extraction cores.}
\item{derivativeStep}{Positive finite relative central-difference step.}
}
\description{Reconstructs study DRIFT modes and their joint delta-method covariance from the optimizer estimate and full raw covariance.}
\details{Native model evaluation preserves nonlinear transformations. These estimates and SEs are distinct from means and SDs of approximate parameter draws. Profiles follow their first occurrence in the fitted data; studyDesign identifies them. Additional moderators and Bayesian fits are not supported. DRIFT estimates are on the original time scale.}
\value{A list with studyEstimates, covariance, jacobian, studyDesign, representativeSubjects (internal subject indices), scaleTime and estimator.}
\seealso{\code{\link{ctmaReMoments}}}
