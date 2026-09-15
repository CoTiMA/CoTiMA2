\name{ctmaFitConfig}
\alias{ctmaFitConfig}
\alias{print.ctmaFitConfig}
\title{Benchmark a CoTiMA fitting configuration}
\usage{
ctmaFitConfig(
  ctmaInitFit = NULL,
  ...,
  backend = c("julia", "stan"),
  backendArgs = list(),
  time_limit = 600,
  candidate_subject_threads = NULL,
  candidate_chunksizes = NULL,
  candidate_carefulFit = c(FALSE, TRUE),
  pilot_iterations = 5L,
  benchmark_repetitions = 2L,
  warmup = TRUE,
  verbose = TRUE
)

\method{print}{ctmaFitConfig}(x, ...)
}
\arguments{
\item{ctmaInitFit}{Object accepted by \code{ctmaFit(ctmaInitFit=...)}.}
\item{...}{Model-building arguments accepted by \code{ctmaFit()}. The
configuration function always uses \code{fit=FALSE}.}
\item{backend}{Estimation backend, \code{"julia"} or \code{"stan"}.}
\item{backendArgs}{Named backend controls. Explicit values are retained as
constraints and recommendations are merged into this list. Calling
\code{ctmaBackendDefaults()} without arguments returns Stan defaults only.}
\item{time_limit}{Approximate pilot budget in seconds. An individual backend
evaluation cannot be interrupted midway.}
\item{candidate_subject_threads}{Optional positive integer Julia concurrency
candidates.}
\item{candidate_chunksizes}{Optional positive ForwardDiff chunk candidates.}
\item{candidate_carefulFit}{Logical values to compare in short Julia optimizer
pilots. By default both joint optimization (\code{FALSE}) and staged
optimization (\code{TRUE}) are tested when the model has time-invariant
effects.}
\item{pilot_iterations}{Positive number of optimizer iterations per stage and
per \code{carefulFit} candidate. These pilots identify a promising strategy;
they are not intended to converge.}
\item{benchmark_repetitions}{Number of repeated likelihood and gradient
evaluations per candidate. Median timings reduce selection by measurement
noise.}
\item{warmup}{Whether to perform one untimed Julia gradient evaluation before
the measured grid.}
\item{verbose}{Whether to print progress and the recommendation.}
\item{x}{A \code{ctmaFitConfig} object.}
}
\value{
An object of class \code{ctmaFitConfig} containing the matching setup,
benchmark table, resolved backend controls, and recommended
\code{backendArgs}.
}
\description{
Creates the same unfitted model setup as \code{ctmaFit()} and performs a
bounded backend pilot. Julia ML pilots rank combinations of runtime-adjustable
subject concurrency and ForwardDiff chunk size at a common raw starting
vector. Unsafe parallel components are excluded from thread comparisons.
Short optimizer pilots compare \code{carefulFit=FALSE} and
\code{carefulFit=TRUE} using both objective improvement per elapsed second and
raw-gradient reduction. For Stan,
supplied controls are validated and retained and setup cost is reported;
pooled ctsem/Stan optimization does not expose the equivalent subject/chunk
tuning grid.
}
\examples{
\dontrun{
config <- ctmaFitConfig(
  ctmaInitFit = initFit,
  backend = "julia",
  time_limit = 600
)
config$recommended
}
}
