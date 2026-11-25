\name{ctmaFitToModel}
\alias{ctmaFitToModel}
\title{Convert a CoTiMA Fit into a Reusable Model}
\usage{ctmaFitToModel(ctmaFitObject, model = NULL, reUseEmprawData = TRUE,
  emprawLong = FALSE)}
\arguments{
\item{ctmaFitObject}{A complete initialization or pooled CoTiMA fit, or a ctmaPrep study list.}
\item{model}{NULL retains the model. Case-insensitive init/ctmaInit and
fit/ctmaFit select standard models. prep/ctmaPrep delegates to ctmaFitToPrep.}
\item{reUseEmprawData}{Reuse cached empirical data. FALSE permits regeneration
from covariance matrices for separate-study models.}
\item{emprawLong}{Additionally include long-format cached data in original units.}
}
\description{Converts fitted study metadata, empirical data and model specifications
into input for ctmaInit or ctmaFit. Stored long times are divided by the saved
scaleTime before it is reapplied exactly once at refitting.}
\details{Explicit fitting arguments override inherited arguments, including
explicit NULL and FALSE. Standard-model conversion discards original constraints,
moderator specifications and study-level random effects. A standard initialization
model disables person-level random intercepts when two-wave studies are present;
explicit indVarying can restore them. Original-model conversion preserves the
separate-study or pooled structure. Convert explicitly to init or fit to change
structure. Missing empirical caches cause an error when reuse is requested.
Prepared-study input generates empirical data once using ctmaInit(fit=FALSE), without fitting primaries. Set the R seed before conversion for reproducible pseudo data. For a preserved pooled model, construction changes other than person variation require conversion to fit.}
\value{A CoTiMAModel, or for prep the result of ctmaFitToPrep.}
\seealso{\code{\link{ctmaFit}}, \code{\link{ctmaInit}}, \code{\link{ctmaFitToPrep}}}
