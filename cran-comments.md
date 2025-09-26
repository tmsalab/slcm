## R CMD check results

0 errors | 0 warnings | 1 note

This release is being made as it relates to changes in RcppArmadillo requiring
packages move to a higher C++ standard.

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘edmdata’
    All declared Imports should be used.

This CRAN check is being extra sensitive regarding a data package that is being
used in both our vignette and examples. The package is listed as an import in the
DESCRIPTION file and is used in the vignette and examples, but not in the main
functions of the package.

* Possibly misspelled words in DESCRIPTION:
  Culpepper (18:28)
  Liang (18:50)
  - Both names correct as they are the author lastnames on the paper describing the SLCM
