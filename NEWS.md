# News

## Version Updates

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
