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
| Microsoft Windows: | [MRMCaov\_0.1.8.zip](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/Ed0_WcuPbe5LkhapFOzmQYMB72GSKzNQmhnpUDWPh7nheg)    |
| Linux:             | [MRMCaov\_0.1.8.tar.gz](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/ESAtW6VQ2MlKu2Dk_dmYbK8B8niyS2p0igKuQ8vtfvqDVA) |

## Citing the Software

``` r
## Text format
citation("MRMCaov")
```

``` 

To cite MRMCaov in publications, please use the following two
references, including the R package URL.

Smith BJ, Hillis SL, Pesce LL (2020). _MCMCaov: Multi-Reader Multi-Case
Analysis of Variance_. R package version 0.1.8, <URL:
https://github.com/brian-j-smith/MRMCaov>.

Smith BJ, Hillis SL (2020). "Multi-reader multi-case analysis of
variance software for diagnostic performance comparison of imaging
modalities." In Samuelson F, Taylor-Phillips S (eds.), _Proceedings of
SPIE 11316, Medical Imaging 2020: Image Perception, Observer
Performance, and Technology Assessment_, 113160K. doi:
10.117/12.2549075 (URL: https://doi.org/10.117/12.2549075), <URL:
https://pubmed.ncbi.nlm.nih.gov/32351258>.

To see these entries in BibTeX format, use 'print(<citation>,
bibtex=TRUE)', 'toBibtex(.)', or set
'options(citation.bibtex.max=999)'.
```

``` r
## Bibtex format
toBibtex(citation("MRMCaov"))
```

    @Manual{MRMCaov-package,
      author = {Brian J Smith and Stephen L Hillis and Lorenzo L Pesce},
      title = {MCMCaov: Multi-Reader Multi-Case Analysis of Variance},
      year = {2020},
      note = {R package version 0.1.8},
      url = {https://github.com/brian-j-smith/MRMCaov},
    }
    
    @InProceedings{MRMCaov-SPIE2020,
      author = {Brian J. Smith and Stephen L. Hillis},
      title = {Multi-reader multi-case analysis of variance software for diagnostic performance comparison of imaging modalities},
      booktitle = {Proceedings of SPIE 11316, Medical Imaging 2020: Image Perception, Observer Performance, and Technology Assessment},
      editor = {Frank Samuelson and Sian Taylor-Phillips},
      month = {16 March},
      year = {2020},
      pages = {113160K},
      doi = {10.117/12.2549075},
      url = {https://pubmed.ncbi.nlm.nih.gov/32351258},
    }

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
