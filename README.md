MRMCaov: Multi-Reader Multi-Case Analysis of Variance
================

## Description

**MRMCaov** is an R package for statistical comparison of diagnostic
tests - such as those based on medical imaging - for which ratings have
been obtained from multiple readers and on multiple cases. Features of
the package include the following.

  - Comparison of imaging modalities or diagnostic tests with respect to
    ROC metrics
  - ROC metrics include AUC, partial AUC, sensitivity, and specificity
    as well as user-defined metrics
  - Parametric and nonparametric estimation and plotting of ROC curves
  - Support for factorial and nested study designs
  - DeLong, jackknife, and unbiased covariance estimation
  - Inference for random or fixed readers and cases
  - Compatibility with Microsoft Windows, Mac OS, and Linux

## Download and Installation

To install the package, do the following:

1.  Download the archive file for your operating system from the links
    provided in the table below.
2.  Run [RStudio](https://www.rstudio.com/products/rstudio/).
3.  Click on “Tools” from the RStudio menu and then select “Install
    Packages…” to bring up the “Install Packages” dialog box.
4.  From the dialog box, select “Package Archive File (.zip; .tar.gz)”
    as the source from which to install, and browse to and select the
    downloaded package archive file.
5.  Install packages `ggplot2`, `gtools`, and `pROC` from the CRAN
    selection in the dialog box if you have not installed these packages
    previously. Alternatively, these can be installed from the R console
    with the command.

<!-- end list -->

``` r
install.packages(c("ggplot2", "gtools", "pROC"))
```

| Operating System   | Package Archive File                                                                                                                    |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------- |
| Microsoft Windows: | [MRMCaov\_0.1.5.zip](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/EXJg0rqU_m5MuXcOB4rpYo0Ba1XZ_J-5JQgr307u-c6EPA)    |
| Non-Windows:       | [MRMCaov\_0.1.5.tar.gz](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/EctXK89M5DRPo3v6r4poCQkB3Y5dUndP3c9ugZxw6_tCEg) |

## Getting Started

``` r
## Load the package
library(MRMCaov)

## View main help page
?mrmc

## Analyze the VanDyke dataset
est <- mrmc(empirical_auc(truth, rating), treatment, reader, case, data = VanDyke)
plot(est)
```

![](README_files/figure-gfm/mrmc-1.png)<!-- -->

``` r
summary(est)
```

    ## Factor types: Random Readers and Random Cases
    ## 
    ## Experimental design: factorial 
    ## 
    ## Obuchowski-Rockette variance component and covariance estimates:
    ## 
    ##                      Estimate Correlation
    ## reader           0.0015349993          NA
    ## treatment:reader 0.0002004025          NA
    ## Error            0.0008022883          NA
    ## Cov1             0.0003466137   0.4320314
    ## Cov2             0.0003440748   0.4288668
    ## Cov3             0.0002390284   0.2979333
    ## 
    ## 
    ## ANOVA global test of equal treatment empirical_auc:
    ## 
    ##         MS(T)      MS(T:R)         Cov2         Cov3 Denominator        F df1
    ## 1 0.004796171 0.0005510306 0.0003440748 0.0002390284 0.001076263 4.456319   1
    ##        df2    p-value
    ## 1 15.25967 0.05166569
    ## 
    ## 
    ## 95% CIs and tests for treatment empirical_auc pairwise differences:
    ## 
    ##   Comparison    Estimate     StdErr       df      CI.Lower      CI.Upper
    ## 1      1 - 2 -0.04380032 0.02074862 15.25967 -0.0879594986  0.0003588544
    ##           t    p-value
    ## 1 -2.110999 0.05166569
    ## 
    ## 
    ## 95% treatment empirical_auc CIs (each analysis based only on data for the specified treatment):
    ## 
    ##    Estimate       MS(R)         Cov2     StdErr       df  CI.Lower  CI.Upper
    ## 1 0.8970370 0.003082629 0.0004839618 0.03317360 12.74465 0.8252236 0.9688505
    ## 2 0.9408374 0.001304602 0.0002041879 0.02156637 12.71019 0.8941378 0.9875369
