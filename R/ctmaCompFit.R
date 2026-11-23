#' ctmaCompFit
#'
#' @description Performs log-liklihood ratio tests to compare the fit of 2 models (CoTiMAFit objects created with \code{\link{ctmaFit}}
#' or \code{\link{ctmaEqual}}), i.e., the difference between the two -2 times LLs between the first model and the more constrained second
#' model. The nested structure of the two models is assumed to be given and not checked.
#'
#' @details Retained numeric fit statistics are used even when a native pooled fit has no legacy model label. Negative minus-two log likelihoods are compared arithmetically. The same observations, nesting, convergence and regular identification assumptions are not verified by this function.
#'
#' @param model1 Model 1
#' @param model2 Model 2
#'
#' @importFrom stats pchisq
#'
#' @examples
#' ctmaLoadExamples(c("CoTiMAFullInv23Fit_6", "CoTiMAFullInvEq23Fit_6"))
#' minus2llDiffTest <- ctmaCompFit(CoTiMAFullInv23Fit_6,
#'                                 CoTiMAFullInvEq23Fit_6)
#' summary(minus2llDiffTest)
#' ctmaUnloadExamples(c("CoTiMAFullInv23Fit_6", "CoTiMAFullInvEq23Fit_6"))
#'
#' @export ctmaCompFit
#'
#' @return Returns the the difference between the two -2 times LLs (Diff_Minus2LL), the associated difference in degrees of
#' freedom (Diff_df (= Diff_n.params)), and the probability (prob).
#'
#'
ctmaCompFit <- function(model1=NULL, model2=NULL) {
  first <- .ctmaComparisonStats(model1, "Model 1")
  second <- .ctmaComparisonStats(model2, "Model 2")
  n.model1 <- length(first$minus2ll)
  n.model2 <- length(second$minus2ll)
  names.model1 <- first$model
  names.model2 <- second$model
  allResults <- list()
  counter <- 0
  for (i in seq_len(n.model1)) {
    for (j in seq_len(n.model2)) {
      counter <- counter + 1; counter
      allResults[[counter]] <- list()
      if (counter > 1) {
        allResults[[counter]]$" " <- "############################    MODEL COMPARISON   ############################"
      } else {
          allResults[[counter]]$" " <- "##########################   NEXT MODEL COMPARISON   ##########################"
        }
      allResults[[counter]]$"Model: " <- paste0(names.model1[i], " COMPARED WITH   Model: ", names.model2[j]); allResults[[counter]]$comparison
      allResults[[counter]]$"Diff_Minus2LL: " <- abs(first$minus2ll[i] - second$minus2ll[j]); allResults[[counter]]$Diff_Minus2LL
      allResults[[counter]]$"Diff_df (= Diff_n.params): " <- abs(first$n.parameters[i] - second$n.parameters[j]); allResults[[counter]]$Diff_df
      if (allResults[[counter]]$"Diff_df (= Diff_n.params): " > 0 ) {
        allResults[[counter]]$"prob: " <- 1-stats::pchisq(allResults[[counter]]$Diff_Minus2LL, allResults[[counter]]$Diff_df); allResults[[counter]]$prob
      } else {
        allResults[[counter]]$"prob: " <- "Diff_df = 0, thus probability cannot be computed"
      }
      allResults[[counter]]$"Message1: " <- "A prob value < .05 indicates a significant difference."; allResults[[counter]]$Message
      allResults[[counter]]$"Message2: " <- "This is valid only if models are nested, which is NOT checked"; allResults[[counter]]$Message
    }
  }
  return(list(paste(names(unlist(allResults)), unlist(allResults))))
}


# Native pooled containers do not have the legacy summary$model label. Extract
# retained fit statistics, never parse the abbreviated print method.
.ctmaComparisonStats <- function(object, label) {
  if (!is.list(object)) stop("Supply two fitted CoTiMA objects.", call.=FALSE)
  s <- object$summary
  if (!is.list(s)) s <- list()
  ll <- s$minus2ll
  k <- s$n.parameters
  if (!length(ll) && is.list(object$estimate)) ll <- -2*object$estimate$loglik
  if (!length(k) && is.list(object$estimate)) k <- length(object$estimate$raw)
  ll <- as.numeric(unlist(ll)); k <- as.numeric(unlist(k))
  if (!length(ll) || length(k)!=length(ll) || any(!is.finite(ll)) ||
      any(!is.finite(k)) || any(k<0)) stop(
    "Each fit must supply matching finite -2 log likelihoods and parameter counts.", call.=FALSE)
  model <- as.character(s$model)
  if (!length(model)) model <- if(length(ll)==1L) label else paste(label,seq_along(ll))
  if(length(model)!=length(ll)) stop("Model labels do not match the fit statistics.",call.=FALSE)
  list(model=model, minus2ll=ll,n.parameters=k)
}
