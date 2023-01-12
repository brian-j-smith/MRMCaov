MRMCaov: Multi-Reader Multi-Case Analysis of Variance
================

[![Generic
badge](https://img.shields.io/badge/docs-online-green.svg)](https://brian-j-smith.github.io/MRMCaov/)

# Getting Started

**MRMCaov** is an R package for statistical comparison of diagnostic
tests - such as those based on medical imaging - for which ratings have
been obtained from multiple readers and on multiple cases. Features of
the package include the following.

- Statistical comparisons of diagnostic tests with respect to reader
  performance metrics
- Comparisons based on the ANOVA model of Obuchowski-Rockette and the
  unified framework of Hillis
- Reader performance metrics for area under receiver operating
  characteristic curves (ROC AUCs), partial AUCs, expected utility of
  ROC curves, likelihood ratio of positive or negative tests,
  sensitivity, specificity, and user-defined metrics
- Parametric and nonparametric estimation and plotting of ROC curves
- Support for factorial, nested, and partially paired study designs
- Inference for random or fixed readers and cases
- Conversion of Obuchowski-Rockette to Roe, Metz & Hillis model
  parameters and vice versa
- DeLong, jackknife, and unbiased covariance estimation
- Compatibility with Microsoft Windows, MacOS, and Linux

## Documentation: [User Guide](https://brian-j-smith.github.io/MRMCaov/using.html)

## Installation

Enter the following at the R console to install the package.

``` r
install.packages("MRMCaov")
```

## Citing the Software

``` r
## Text format
citation("MRMCaov")
```


    To cite MRMCaov in publications, please use the following two references, including the R package URL.

      Smith BJ, Hillis SL, Pesce LL (2023). _MCMCaov: Multi-Reader Multi-Case Analysis of Variance_. R
      package version 0.3.0, <https://github.com/brian-j-smith/MRMCaov>.

      Smith BJ, Hillis SL (2020). "Multi-reader multi-case analysis of variance software for diagnostic
      performance comparison of imaging modalities." In Samuelson F, Taylor-Phillips S (eds.),
      _Proceedings of SPIE 11316, Medical Imaging 2020: Image Perception, Observer Performance, and
      Technology Assessment_, 113160K. doi:10.1117/12.2549075 <https://doi.org/10.1117/12.2549075>,
      <https://pubmed.ncbi.nlm.nih.gov/32351258>.

    To see these entries in BibTeX format, use 'print(<citation>, bibtex=TRUE)', 'toBibtex(.)', or set
    'options(citation.bibtex.max=999)'.

``` r
## Bibtex format
toBibtex(citation("MRMCaov"))
```

    @Manual{MRMCaov-package,
      author = {Brian J Smith and Stephen L Hillis and Lorenzo L Pesce},
      title = {{MCMCaov}: Multi-Reader Multi-Case Analysis of Variance},
      year = {2023},
      note = {R package version 0.3.0},
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
      doi = {10.1117/12.2549075},
      url = {https://pubmed.ncbi.nlm.nih.gov/32351258},
    }
