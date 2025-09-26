

<!-- README.md is generated from README.qmd. Please edit that file -->

# slcm

<!-- badges: start -->

[![R-CMD-check](https://github.com/tmsalab/slcm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tmsalab/slcm/actions/workflows/R-CMD-check.yaml)
[![Package-License](http://img.shields.io/badge/license-GPL%20(%3E=2)-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html)
<!-- badges: end -->

The goal of `slcm` is to provide an implementation of the exploratory
Sparse Latent Class Model (SLCM) for Binary Data described by Chen, Y.,
Culpepper, S. A., and Liang, F. (2020) <doi:10.1007/s11336-019-09693-2>.

This package contains a new implementation of the proposed SLCM based on
the paper. You may find original papers implementation in the [`inst/`
folder](https://github.com/tmsalab/slcm/tree/main/inst) of the package.

## Installation

You can install the released version of slcm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("slcm")
```

Or, you can be on the cutting-edge development version on
[GitHub](https://github.com/) using:

``` r
# install.packages("devtools")
devtools::install_github("tmsalab/slcm")
```

## Usage

To use `slcm`, load the package using:

``` r
library("slcm")
```

From here, the SLCM model can be estimated using:

``` r
model_slcm = slcm::slcm(
  y = <data>,
  k = <k>
)
```

## Authors

James Joseph Balamuta and Steven Andrew Culpepper

## Citing the `slcm` package

To ensure future development of the package, please cite `slcm` package
if used during an analysis or simulation study. Citation information for
the package may be acquired by using in *R*:

``` r
citation("slcm")
```

## License

GPL (\>= 2)
