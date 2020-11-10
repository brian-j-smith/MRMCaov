MRMCaov: Multi-Reader Multi-Case Analysis of Variance
================

[![Generic
badge](https://img.shields.io/badge/docs-online-green.svg)](https://brian-j-smith.github.io/MRMCaov/)

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
5.  Install the packages below with the given command submitted at the R
    console, if they have not been installed previously.

<!-- end list -->

``` r
install.packages(c("ggplot2", "gtools", "mvtnorm", "pROC", "tibble"))
```

| Operating System   | Package Archive File                                                                                                                     |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------- |
| Microsoft Windows: | [MRMCaov\_0.1.12.zip](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/EZvSfeD6QEBKstoYf3ffWDwBlPS_CNBhuTkqKEGcHZt93w)    |
| Linux:             | [MRMCaov\_0.1.12.tar.gz](https://iowa-my.sharepoint.com/:u:/g/personal/bjsmith_uiowa_edu/ERjzWoVL-E9HkW893DLwY3EBD-PQzGK8iBQ9yzN7ocHt2Q) |

## Citing the Software

``` r
## Text format
citation("MRMCaov")
```

``` 

To cite MRMCaov in publications, please use the following two references, including the R
package URL.

Smith BJ, Hillis SL, Pesce LL (2020). _MCMCaov: Multi-Reader Multi-Case Analysis of
Variance_. R package version 0.1.12, <URL: https://github.com/brian-j-smith/MRMCaov>.

Smith BJ, Hillis SL (2020). "Multi-reader multi-case analysis of variance software for
diagnostic performance comparison of imaging modalities." In Samuelson F, Taylor-Phillips
S (eds.), _Proceedings of SPIE 11316, Medical Imaging 2020: Image Perception, Observer
Performance, and Technology Assessment_, 113160K. doi: 10.117/12.2549075 (URL:
https://doi.org/10.117/12.2549075), <URL: https://pubmed.ncbi.nlm.nih.gov/32351258>.

To see these entries in BibTeX format, use 'print(<citation>, bibtex=TRUE)',
'toBibtex(.)', or set 'options(citation.bibtex.max=999)'.
```

``` r
## Bibtex format
toBibtex(citation("MRMCaov"))
```

    @Manual{MRMCaov-package,
      author = {Brian J Smith and Stephen L Hillis and Lorenzo L Pesce},
      title = {MCMCaov: Multi-Reader Multi-Case Analysis of Variance},
      year = {2020},
      note = {R package version 0.1.12},
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
