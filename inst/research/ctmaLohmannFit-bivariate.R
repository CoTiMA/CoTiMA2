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

# TWO VARIABLES: literal original-supplement model, including its caveats.
# Set this to your local supplement root or its modelcompiler.R file.
supplementDirectory <- "/Users/cdormann/SynologyDrive/Synology Drive/CHRISTIAN/TEXTE/LITERATUR/ARTIKEL/LOHMANN ET AL SUPP MAT"

# Generate valid stationary reference correlations in X1,Y1,X2,Y2,... order.
# These follow the intended model, not the supplement's residual-formula discrepancy.
# Diagonal drift keeps this example simple; estimation still frees cross effects.
R_true <- lapply(seq_along(N), function(s) {
  A <- diag(c(-0.4-0.03*s, -0.6-0.02*s))
  tau <- c(0.2+0.02*s, 0.3+0.02*s)
  T <- diag(tau)
  T[1,2] <- T[2,1] <- 0.2*sqrt(prod(tau))
  Phi <- diag(1-tau)
  Phi[1,2] <- Phi[2,1] <- 0.1*sqrt(prod(1-tau))
  R <- matrix(0, 2*waves, 2*waves)
  for (i in seq_len(waves)) for (j in seq_len(waves)) {
    dt <- abs(i-j)*lags[s]
    E <- diag(exp(diag(A)*dt))
    block <- if (i >= j) T + E %*% Phi else T + Phi %*% E
    R[(2*i-1):(2*i),(2*j-1):(2*j)] <- block
  }
  R
})
R_list <- lapply(seq_along(N), function(s)
  cov2cor(stats::rWishart(1, N[s]-1, R_true[[s]])[,,1]))
V_list <- Map(corSamplingCov, R_list, as.list(N))

specification <- CoTiMA2:::ctmaLohmannFit(
  n.latent = 2,
  correlationMatrices = R_list,
  sampleSizes = N,
  timeLag = lags,
  samplingCovariances = V_list,
  supplementDirectory = supplementDirectory,
  reStructure = "URE", # Also supports FE, ERE and CRE.
  fit = FALSE
)
cat(specification$modelCode)

if (runSampling) {
  fit_bivariate <- CoTiMA2:::ctmaLohmannFit(
    n.latent = 2, stanData = specification$stanData,
    supplementDirectory = supplementDirectory,
    reStructure = "URE", chains = 4, iter = 2000, warmup = 1000,
    coresToUse = 4, seed = 20261106,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  print(fit_bivariate$summary)
  rstan::check_hmc_diagnostics(fit_bivariate$fit)
}
# Five studies illustrate syntax, not adequate RE precision or MC recovery.
