\name{pooledJuliaFit}
\alias{pooledJuliaFit}
\alias{print.CoTiMAPooledJuliaFit}
\alias{summary.CoTiMAPooledJuliaFit}
\alias{coef.CoTiMAPooledJuliaFit}
\alias{vcov.CoTiMAPooledJuliaFit}
\title{Methods for pooled Julia CoTiMA fits}
\usage{
\method{print}{CoTiMAPooledJuliaFit}(x, ...)
\method{summary}{CoTiMAPooledJuliaFit}(object, ...)
\method{coef}{CoTiMAPooledJuliaFit}(object, ...)
\method{vcov}{CoTiMAPooledJuliaFit}(object, ...)
}
\arguments{
\item{x, object}{A pooled Julia CoTiMA fit.}
\item{...}{Arguments forwarded to the corresponding native-fit method.}
}
\description{Pooled Julia results returned by ctmaFit() retain the full native fit in studyFitList. Summaries retain the native fields and additionally report
transformed DRIFT and DIFFUSIONcov estimates on the original time scale.
Raw parameters and diagnostics are not rescaled. Other methods delegate
to the unchanged engine object.}
\value{The native summary plus \code{estimates_fitted_time_scale},
\code{estimates_original_time_scale} and \code{time_scale_info}; or the native
coefficients or covariance matrix. Native DT intervals retain their fitted
time units; \code{time_scale_info} records the corresponding original intervals. Printing returns the pooled object invisibly.}
\seealso{\code{\link{ctmaFit}}, \code{\link{ctmaFitDiagnostics}}}
