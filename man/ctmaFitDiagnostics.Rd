\name{ctmaFitDiagnostics}
\alias{ctmaFitDiagnostics}
\title{Separate convergence, identification and scale diagnostics}
\usage{ctmaFitDiagnostics(object, print = TRUE)}
\arguments{
\item{object}{A CoTiMA2 Julia fit or a conventional CoTiMAFit object.}
\item{print}{Whether to print the diagnostic report.}}
\description{Cached Laplace integration-accuracy status is reported separately from optimizer convergence. Unresolved or missing verdicts do not imply an accuracy pass. Reports read stored metadata without rerunning a check.
Reports stored optimizer, identification, transform-saturation and
posterior diagnostics without refitting or starting Julia. Missing fields are unknown.}
\details{New native fits cache population matrices at fitting time. Binary
linear-predictor SD above 10 at Time 0 or at the population stationary covariance
is a labelled heuristic, not an identification test, validated numerical cutoff
or instruction to fix diffusion. No parameter or standard error is changed.
Optimizer convergence does not certify identification, reliable standard errors,
likelihood accuracy or posterior convergence. Older saved fits can be examined
for their stored diagnostics; raw-coordinate and retained transformed-summary
diagnostics are labelled separately. SMC R-hat is not applicable. Native pooled
fits and primary-study rosters are recognized independently; unavailable matrix checks are reported as not assessed.}
\value{Invisibly, a diagnostic list, or study-indexed lists for ctmaInit fits.}
