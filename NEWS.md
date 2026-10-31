# slcm 0.1.1.9000

## Bug fixes

- `slcm()` now stores every draw when `burnin = 0`. Previously the
  comparison `t > burnin - 1` wrapped around on an unsigned integer, so
  nothing was stored and the returned chains and the `delta` and `q`
  estimates were all zeros.
- `slcm()` now stops with an error when `y` contains missing values.
  Previously a single `NA` produced an invalid fit without warning.

## Internal

- Added a testthat test suite.
