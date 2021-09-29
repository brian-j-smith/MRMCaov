#' @details
#' The following functions are available in \pkg{MRMCaov} for estimation and
#' comparison of test performance metrics in studies involving multiple cases
#' and one or more readers.
#'
#' Statistical Inference:
#' \tabular{ll}{
#'   \code{\link{mrmc}} \tab Multi-reader multi-case ANOVA \cr
#'   \code{\link{srmc}} \tab Single-reader multi-case ANOVA \cr
#'   \code{\link{stmc}} \tab Single-test (single-reader) multi-case Estimation \cr
#' }
#'
#' Tabular and Graphical Summaries:
#' \tabular{ll}{
#'   \code{\link{parameters}} \tab ROC curve parameters \cr
#'   \code{\link{plot}}       \tab ROC curve plots \cr
#'   \code{\link{roc_curves}} \tab ROC curves \cr
#'   \code{\link{summary}}    \tab Statistical analysis summaries \cr
#' }
#'
#' Performance Metrics (Binary Rating):
#' \tabular{ll}{
#'   \code{\link{binary_sens}}      \tab Sensitivity \cr
#'   \code{\link{binary_spec}}      \tab Specificity \cr
#' }
#'
#' Performance Metrics (Ordinal or Numeric Rating):
#' \tabular{ll}{
#'   \code{\link{binormal_auc}}     \tab Binormal ROC AUC \cr
#'   \code{\link{binormal_sens}}    \tab ... sensitivity \cr
#'   \code{\link{binormal_spec}}    \tab ... specificity \cr
#'   \code{\link{binormalLR_auc}}   \tab Binormal likelihood ratio ROC AUC \cr
#'   \code{\link{binormalLR_sens}}  \tab ... sensitivity \cr
#'   \code{\link{binormalLR_spec}}  \tab ... specificity \cr
#'   \code{\link{empirical_auc}}    \tab Empirical ROC AUC \cr
#'   \code{\link{empirical_sens}}   \tab ... sensitivity \cr
#'   \code{\link{empirical_spec}}   \tab ... specificity \cr
#'   \code{\link{trapezoidal_auc}}  \tab Empirical ROC AUC \cr
#'   \code{\link{trapezoidal_sens}} \tab ... sensitivity \cr
#'   \code{\link{trapezoidal_spec}} \tab ... sensitivity \cr
#' }
#'
#' Performance Metric Covariance Estimation Methods:
#' \tabular{l}{
#'   \code{\link{DeLong}} \cr
#'   \code{\link{jackknife}} \cr
#'   \code{\link{unbiased}} \cr
#' }
#'
#' ROC Curves:
#' \tabular{ll}{
#'   \code{\link{roc_curves}} \tab Estimate one or more curves \cr
#'   \code{\link{parameters}} \tab Extract curve parameters \cr
#'   \code{\link{points}}     \tab Extract curve points \cr
#'   \code{\link{mean}}       \tab Compute the mean of multiple curves \cr
#'   \code{\link{plot}}       \tab Plot curves \cr
#' }
#'
#' @note
#' This research was supported by the National Institute of Biomedical Imaging
#' and Bioengineering (NIBIB) of the National Institutes of Health under Award
#' Number R01EB025174
#'
"_PACKAGE"
