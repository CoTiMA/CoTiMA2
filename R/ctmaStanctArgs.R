#' Legacy Stan backend defaults
#'
#' @name CoTiMAStanctArgs
#' @export CoTiMAStanctArgs
#' @description Deprecated compatibility snapshot of the Stan backend
#' defaults. New code should use \code{ctmaBackendDefaults("stan")} and pass
#' modifications through \code{backendArgs}.
#'
.ctmaStanBackendDefaults <- function() {
  list(test=TRUE,
                       scaleTI=TRUE, scaleMod=FALSE, scaleClus=FALSE,
                       scaleTD=FALSE, scaleLongData=FALSE, scaleTime=1/1,
                       scaleTimeAuto=NULL,
                       datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                       iter=1000, intoverstates=TRUE,
                       binomial=FALSE,
                       fit=TRUE,
                       intoverpop="auto",
                       sameInitialTimes=FALSE,
                       stationary=FALSE,
                       plot=FALSE, derrind=NA,
                       optimize=TRUE,
                       optimcontrol=list(init = "random",
                                         initsd = 0.01,
                                         #sampleinit = NA,
                                         #deoptim = FALSE,
                                         estonly = FALSE,
                                         tol = 1e-08,
                                         #decontrol = list(),
                                         stochastic = TRUE,
                                         priors = TRUE,
                                         carefulfit = TRUE,
                                         #bootstrapUncertainty = FALSE,
                                         subsamplesize = 1,
                                         parsteps = c(),
                                         parstepsAutoModel = FALSE,
                                         groupFreeThreshold = 0.5,
                                         plot = FALSE,
                                         #is = FALSE,
                                         #isloopsize = 1000,
                                         #isitersize = 1000,
                                         #isESS = 100,
                                         finishsamples = 1000,
                                         #tdf = 10,
                                         #chancethreshold = 100,
                                         #finishmultiply = 5,
                                         lproughnesstarget = 0.2,
                                         cores = 2,
                                         matsetup = NA,
                                         nsubsets = 100,
                                         stochasticTolAdjust = 1000),
                       nlcontrol=list(),
                       #nopriors=TRUE,
                       priors=FALSE,
                       chains=2,
                       #resample=NULL,
                       chais=2,
                       cores=2,
                       # NULL means half of `iter`, matching rstan's default.
                       # Supplying backendArgs$warmup overrides this value.
                       warmup=NULL,
                       inits=NULL,
                       forcerecompile=FALSE,
                       saveCompile=TRUE,
                       savescores=FALSE,
                       savesubjectmatrices=FALSE,
                       saveComplexPars = FALSE,
                       gendata=FALSE,
                       control=list(adapt_delta = .8,
                                    adapt_window=25, #adapt_window=2,
                                    max_treedepth=10,
                                    adapt_init_buffer=75, #adapt_init_buffer=2,
                                    stepsize = .001),
                       verbose = 0,
                       vb=FALSE)
}

# Keep a snapshot for source compatibility. Internal code must obtain a fresh
# list from .ctmaStanBackendDefaults(), so user modification of this exported
# object cannot change package defaults.
CoTiMAStanctArgs <- .ctmaStanBackendDefaults()

.ctmaStanWarmup <- function(arguments) {
  iterations <- suppressWarnings(as.integer(arguments$iter)[1L])
  warmup <- if (is.null(arguments$warmup)) {
    iterations %/% 2L
  } else {
    suppressWarnings(as.integer(arguments$warmup)[1L])
  }
  if (is.na(iterations) || iterations < 1L || is.na(warmup) || warmup < 0L ||
      warmup >= iterations) {
    stop(
      "Stan sampling requires a positive 'iter' and 'warmup' in the range ",
      "0 <= warmup < iter.",
      call. = FALSE
    )
  }
  warmup
}
