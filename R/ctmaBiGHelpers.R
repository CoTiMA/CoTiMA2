.ctmaBiGMetaAnalysis <- function(effects, standard_errors,
    confidence = 0.95) {
  effects <- as.matrix(effects)
  standard_errors <- as.matrix(standard_errors)
  if (!identical(dim(effects), dim(standard_errors))) {
    stop("Effects and standard errors must have identical dimensions.",
      call. = FALSE)
  }
  confidence <- as.numeric(confidence)[1L]
  if (!is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("'confidence' must be between zero and one.", call. = FALSE)
  }
  z_critical <- stats::qnorm(1 - (1 - confidence) / 2)
  parameter_names <- colnames(effects)
  if (is.null(parameter_names)) {
    parameter_names <- paste0("parameter_", seq_len(ncol(effects)))
  }

  common_names <- c(
    "MeanOfDriftValues", "FixedEffect_Drift",
    "FixedEffect_DriftVariance", "FixedEffect_DriftSE",
    "FixedEffect_DriftUpperLimit", "FixedEffect_DriftLowerLimit",
    "FixedEffect_DriftZ", "FixedEffect_DriftProb"
  )
  heterogeneity_names <- c(
    "AnalysedStudies", "Heterogeneity_df", "tau2Drift", "Q_Drift",
    "H2_Drift", "H2DriftUpperLimit", "H2DriftLowerLimit", "I2_Drift",
    "I2DriftUpperLimit", "I2DriftLowerLimit"
  )
  random_names <- c(
    "RandomEffecttot_Drift", "RandomEffecttot_DriftVariance",
    "RandomEffecttot_DriftSE", "RandomEffecttot_DriftUpperLimit",
    "RandomEffecttot_DriftLowerLimit", "RandomEffecttot_DriftZ",
    "RandomEffecttot_DriftProb", "RandomEffecttot_DriftUpperLimitPI",
    "RandomEffecttot_DriftLowerLimitPI"
  )

  analyse_one <- function(effect, standard_error) {
    valid <- is.finite(effect) & is.finite(standard_error) &
      standard_error > 0
    effect <- effect[valid]
    standard_error <- standard_error[valid]
    k <- length(effect)
    if (k < 2L) {
      common <- stats::setNames(rep(NA_real_, length(common_names)),
        common_names)
      heterogeneity <- stats::setNames(
        rep(NA_real_, length(heterogeneity_names)), heterogeneity_names)
      heterogeneity["AnalysedStudies"] <- k
      random <- stats::setNames(rep(NA_real_, length(random_names)),
        random_names)
      return(list(k = k, common = common,
        heterogeneity = heterogeneity, random = random))
    }

    weights <- 1 / standard_error^2
    sum_weights <- sum(weights)
    common_mean <- sum(weights * effect) / sum_weights
    common_variance <- 1 / sum_weights
    common_se <- sqrt(common_variance)
    common_z <- common_mean / common_se
    common_p <- 2 * stats::pnorm(-abs(common_z))
    q <- sum(weights * (effect - common_mean)^2)
    df <- k - 1L
    h2 <- max(1, q / df)
    i2 <- max(0, (q - df) / q) * 100
    if (!is.finite(i2)) i2 <- 0

    c_value <- sum_weights - sum(weights^2) / sum_weights
    tau2 <- if (is.finite(c_value) && c_value > 0) {
      max(0, (q - df) / c_value)
    } else NA_real_

    se_log_h <- if (k > 2L) {
      if (q > k) {
        0.5 * (log(q) - log(df)) /
          (sqrt(2 * q) - sqrt(2 * k - 3))
      } else {
        sqrt((1 / (2 * (k - 2))) *
          (1 - 1 / (3 * (k - 2)^2)))
      }
    } else NA_real_
    h <- sqrt(h2)
    h_lower <- if (is.finite(se_log_h)) {
      exp(log(h) - z_critical * se_log_h)
    } else NA_real_
    h_upper <- if (is.finite(se_log_h)) {
      exp(log(h) + z_critical * se_log_h)
    } else NA_real_
    h2_lower <- h_lower^2
    h2_upper <- h_upper^2
    i2_lower <- if (is.finite(h_lower)) max(0, (h_lower^2 - 1) / h_lower^2 * 100) else NA_real_
    i2_upper <- if (is.finite(h_upper)) max(0, (h_upper^2 - 1) / h_upper^2 * 100) else NA_real_

    if (is.finite(tau2)) {
      random_weights <- 1 / (standard_error^2 + tau2)
      random_mean <- sum(random_weights * effect) / sum(random_weights)
      random_variance <- 1 / sum(random_weights)
      random_se <- sqrt(random_variance)
      random_z <- random_mean / random_se
      random_p <- 2 * stats::pnorm(-abs(random_z))
      prediction_critical <- if (k > 2L) {
        stats::qt(1 - (1 - confidence) / 2, df = k - 2L)
      } else NA_real_
      prediction_half_width <- prediction_critical *
        sqrt(tau2 + random_variance)
    } else {
      random_mean <- random_variance <- random_se <- random_z <-
        random_p <- prediction_half_width <- NA_real_
    }

    list(
      k = k,
      common = c(
        MeanOfDriftValues = mean(effect),
        FixedEffect_Drift = common_mean,
        FixedEffect_DriftVariance = common_variance,
        FixedEffect_DriftSE = common_se,
        FixedEffect_DriftUpperLimit = common_mean + z_critical * common_se,
        FixedEffect_DriftLowerLimit = common_mean - z_critical * common_se,
        FixedEffect_DriftZ = common_z,
        FixedEffect_DriftProb = common_p
      ),
      heterogeneity = c(
        AnalysedStudies = k,
        Heterogeneity_df = df,
        tau2Drift = tau2,
        Q_Drift = q,
        H2_Drift = h2,
        H2DriftUpperLimit = h2_upper,
        H2DriftLowerLimit = h2_lower,
        I2_Drift = i2,
        I2DriftUpperLimit = i2_upper,
        I2DriftLowerLimit = i2_lower
      ),
      random = c(
        RandomEffecttot_Drift = random_mean,
        RandomEffecttot_DriftVariance = random_variance,
        RandomEffecttot_DriftSE = random_se,
        RandomEffecttot_DriftUpperLimit = random_mean + z_critical * random_se,
        RandomEffecttot_DriftLowerLimit = random_mean - z_critical * random_se,
        RandomEffecttot_DriftZ = random_z,
        RandomEffecttot_DriftProb = random_p,
        RandomEffecttot_DriftUpperLimitPI = random_mean + prediction_half_width,
        RandomEffecttot_DriftLowerLimitPI = random_mean - prediction_half_width
      )
    )
  }

  analyses <- lapply(seq_len(ncol(effects)), function(column) {
    analyse_one(effects[, column], standard_errors[, column])
  })
  bind_component <- function(component) {
    out <- do.call(cbind, lapply(analyses, `[[`, component))
    colnames(out) <- parameter_names
    out
  }
  list(
    common = bind_component("common"),
    heterogeneity = bind_component("heterogeneity"),
    random = bind_component("random"),
    k = vapply(analyses, `[[`, integer(1L), "k")
  )
}
