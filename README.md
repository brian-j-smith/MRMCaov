MRMCaov: Multi-Reader Multi-Case Analysis of Variance
================

## Description

**MRMCaov** is an R package for statistical comparison of diagnostic
tests - such as those based on medical imaging - for which ratings have
been obtained from multiple readers and on multiple cases. Features of
the package include the following.

  - Performance metrics: area under the receiver operating
    characteristic curve (empirical/trapezoidal and proper parametric
    curves)
  - Study designs: factorial, cases nested within readers, cases nested
    within tests, partially-paired cases
  - Covariance estimation methods: DeLong, jackknife, unbiased

## Download and Installation

To install the package, download the archive file for your operating
system from the links provided in the table below. Next, do the
following:

1.  Run [RStudio](https://www.rstudio.com/products/rstudio/).
2.  Click on “Tools” from the RStudio menu and then select “Install
    Packages…” to bring up the “Install Packages” dialog box.
3.  From the dialog box, select “Package Archive File (.zip; .tar.gz)”
    as the source from which to install, and browse to and select the
    downloaded package archive
file.

| Operating System   | Package Archive File                                                                                                                             |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Microsoft Windows: | [MRMCaov\_0.1.2.zip](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/EZlO1rYXXZ9HiT4fD4Q5AF8B-tyJQw087MlMfbsp-4nZEQ?e=mvIZPk)    |
| Non-Windows:       | [MRMCaov\_0.1.2.tar.gz](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/EVvhUJr9mxJNixFiQVvFx8oB11kISQDxWz_Ay4gIPWb1nw?e=ZmIlUf) |

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
    ##         MS(T)      MS(T:R)         Cov2         Cov3 Denominator        F
    ## 1 0.004796171 0.0005510306 0.0003440748 0.0002390284 0.001076263 4.456319
    ##   df1      df2    p-value
    ## 1   1 15.25967 0.05166569
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
    ##    Estimate       MS(R)         Cov2     StdErr       df  CI.Lower
    ## 1 0.8970370 0.003082629 0.0004839618 0.03317360 12.74465 0.8252236
    ## 2 0.9408374 0.001304602 0.0002041879 0.02156637 12.71019 0.8941378
    ##    CI.Upper
    ## 1 0.9688505
    ## 2 0.9875369
