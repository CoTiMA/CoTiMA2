# CoTiMA2 internal research examples. No sampling unless runSampling is TRUE.
# Run with the newly installed CoTiMA2 package (version 0.14.59 or newer).
# This function is not exported: use CoTiMA2:::ctmaLohmannFit().
runSampling <- FALSE
set.seed(20261106)
N <- rep(500, 5)
lags <- c(0.5, 1, 1.5, 2, 3)
waves <- 4

# Gaussian asymptotic covariance of sample correlations (plug-in estimate).
# Ordering matches R[lower.tri(R)]. This is an explicit choice, not the
# pooled-correlation approximation used in the authors' empirical example.
corSamplingCov <- function(R, n) {
  pairs <- which(lower.tri(R), arr.ind = TRUE)
  V <- matrix(0, nrow(pairs), nrow(pairs))
  for (a in seq_len(nrow(pairs))) for (b in seq_len(nrow(pairs))) {
    i <- pairs[a, 1]; j <- pairs[a, 2]
    k <- pairs[b, 1]; l <- pairs[b, 2]
    V[a,b] <- (0.5*R[i,j]*R[k,l]*(R[i,k]^2+R[i,l]^2+R[j,k]^2+R[j,l]^2) +
      R[i,k]*R[j,l]+R[i,l]*R[j,k] -
      R[i,j]*(R[i,k]*R[i,l]+R[j,k]*R[j,l]) -
      R[k,l]*(R[i,k]*R[j,k]+R[i,l]*R[j,l]))/n
  }
  V
}

# ONE VARIABLE: new univariate extension, not a model from the supplement.
# Var(observed) = trait variance + stationary within-person variance = 1.
# Q = -2*A*(1-trait); Q is determined, not separately free.
rates <- c(-0.35, -0.45, -0.50, -0.55, -0.65)
traits <- c(0.25, 0.30, 0.35, 0.40, 0.45)
R_true <- lapply(seq_along(N), function(s) {
  distance <- abs(outer(0:(waves-1), 0:(waves-1), "-"))*lags[s]
  traits[s] + (1-traits[s])*exp(rates[s]*distance)
})
# Wishart sampling supplies noisy, positive-definite sample correlations.
R_list <- lapply(seq_along(N), function(s)
  cov2cor(stats::rWishart(1, N[s]-1, R_true[[s]])[,,1]))
V_list <- Map(corSamplingCov, R_list, as.list(N))

specification <- CoTiMA2:::ctmaLohmannFit(
  n.latent = 1,
  correlationMatrices = R_list,
  sampleSizes = N,
  timeLag = lags,
  samplingCovariances = V_list,
  reStructure = "URE", # Also supports FE, ERE and CRE.
  fit = FALSE
)
# Inspect the exact specification and prepared data without compiling.
cat(specification$modelCode)

if (runSampling) {
  fit_univariate <- CoTiMA2:::ctmaLohmannFit(
    n.latent = 1, stanData = specification$stanData,
    reStructure = "URE", chains = 4, iter = 2000, warmup = 1000,
    coresToUse = 4, seed = 20261106,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  print(fit_univariate$summary)
  rstan::check_hmc_diagnostics(fit_univariate$fit)
}
# Five studies illustrate syntax, not adequate RE precision or MC recovery.
