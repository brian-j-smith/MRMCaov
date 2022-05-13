# News

## Version Updates

## 0.2.1
* Order tests and readers lexically in results from `mrmc()`.
* Rename `bichisquare` column to `bichisquared` in `binormalLR_curve`.

## 0.2.0
* Implement fast algorithm for unbiased covariance estimation.
* Remove deprecated functions `proproc_roc()`, `proproc_sens()`, and `proproc_spec()`.
* Rename argument `method` to `cov` in `mrmc()`, `srmc()`, and `stmc()`.

## 0.1.16
* Add functions `OR_to_RMH()` and `RMH_to_OR()` for conversion from Obuchowski-Rockette to Roe, Metz & Hillis model parameters and vice versa.

## 0.1.15
* Implement ROC expected utility metrics `binormal_eu()`, `binormalLR_eu()`, and `empirical_eu()`.

## 0.1.14
* Add progress bar for unbiased covariance calculation.
* Truncate confidence intervals only for `sens`, `spec`, and `auc` metrics.
* Fix direction of ratings in `pROC` calculation of empirical ROC AUC to allow values < 0.5.

## 0.1.13
* Add `mrmc()` support for single readers.
* Rename `proproc` functions to `binormalLR`.
* Add `parameters()` method function for `mrmc` output.
* Add `roc_curves()` method function for `srmc` output.
* Add `normalize` argument to `auc()` metrics.
* Redefine `srmc()` function for single-reader (multi-test) multi-case analysis.
* Add `stmc()` function for single-reader (single-test) multi-case analysis.

## 0.1.12
* Fix `proproc_sens()` and `proproc_spec()` to use proproc instead of empirical curves.

## 0.1.11
* Add support for readers nested within tests.
* Add `smmc()` for single reader multiple case analysis.
* Add `parameters()` function to extract ROC curve parameters from `roc_curves()` results.
* Fix partial AUC calculation when specificity ranges are given.
* Add plot option to overlay empirical ROC points on parametric curves.

## 0.1.10
* Subset data in calls to `mrmc()` by event/non-event cases for `binary_sens()`/`binary_spec()` metrics.
* Add `coord_fixed` argument to `plot()` methods.
* Add `mean()` method for `binormal_curves` to average over their ROC curve parameters.

## 0.1.9
* Add `binary_sens()` and `binary_spec()` metrics.

## 0.1.8
* Support sensitivity values in `points()` methods.
* Implement conventional binormal ROC curves and metrics.
* Add package citation.

## 0.1.7
* Add `roc_curves()` function for construction of ROC curves and method functions `points()` for calculating points on the curves, `mean()` for averaging multiple curves, and `plot()` for displaying them graphically.
* Deprecate `roc()` and `proproc()`.

## 0.1.6
* Fix unbiased covariance calculation.

## 0.1.5
* Internal code changes

## 0.1.4
* Internal code changes.

## 0.1.3
* Add metric functions for sensitivity, specificity, and partial AUC.

## 0.1.2
* Initial GitHub upload.
