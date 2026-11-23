\name{ctmaLohmannFit}
\alias{ctmaLohmannFit}
\title{Internal Adapter to the Lohmann et al. Original Stan Models}
\description{Research-only literal replication adapter. Not exported. Reads one original Stan program from a local supplement without executing the compiler script.}
\usage{
ctmaLohmannFit(primaryStudies = NULL, reStructure = "URE", n.latent = 2L,
  supplementDirectory = NULL, correlationMatrices = NULL,
  sampleSizes = NULL, timeLag = NULL, samplingCovariances = NULL,
  stanData = NULL, fit = TRUE, chains = 4L, iter = 2000L,
  warmup = floor(iter / 2), coresToUse = 1L, seed = 20261106L,
  control = list(), init = 0, stanModel = NULL,
  activeDirectory = NULL, ...)
}
\arguments{
\item{primaryStudies}{Optional ctmaPrep study list or convertible complete CoTiMA fit.}
\item{n.latent}{1 for the authored univariate extension; 2 for the literal original bivariate supplement model.}
\item{reStructure}{FE, URE, ERE or CRE, ignoring case.}
\item{supplementDirectory}{Local supplement root or explicit modelcompiler.R file.}
\item{correlationMatrices}{List of complete matrices ordered X1,Y1,X2,Y2,... .}
\item{sampleSizes}{One sample size greater than one per matrix.}
\item{timeLag}{One positive equally spaced interval per study.}
\item{samplingCovariances}{Explicit positive-definite sampling covariance matrices, ordered as each correlation matrix's column-major lower triangle.}
\item{stanData}{Alternatively, the authors' prepared Stan data. Do not combine with study inputs.}
\item{fit}{FALSE prepares code and data without compilation or sampling.}
\item{chains, iter, warmup, coresToUse, seed, control, init}{Stan NUTS settings. Iter includes warmup; at most min(coresToUse, chains) parallel chains.}
\item{stanModel}{Optional compiled rstan model with exactly matching model code.}
\item{activeDirectory}{Optional existing directory for saving the result, overwriting ctmaLohmannFit_DIMD_STRUCTURE.rds.}
\item{...}{Unknown arguments are rejected.}
}
\value{A separate research result with fit, raw rstan summary, Stan data, original model code, source fingerprints and caveats. Not a ctsem fit.}
\details{Access with CoTiMA2:::ctmaLohmannFit(). This deliberately retains the original supplement's residual covariance formula, lack of joint drift/diffusion admissibility checks, and different FE priors/correlation transforms. It is not a certified implementation of the intended stationary mathematical model. No automatic sampling covariance approximation or time rescaling is selected. Full English documentation and examples are installed at inst/research/ctmaLohmannFit.md (system.file("research", "ctmaLohmannFit.md", package="CoTiMA2")).}
\keyword{internal}
