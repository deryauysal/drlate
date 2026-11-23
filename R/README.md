# drlate (R)

R version of the Stata command `drlate`: weighting and doubly robust estimators of the local average treatment effect (LATE) and the local average treatment effect on the treated (LATT), with logit, CBPS or IPT instrument propensity scores.

Reference: Słoczyński, T., Uysal, S. D., and Wooldridge, J. M. "Doubly Robust Estimation of Local Average Treatment Effects Using Inverse Probability Weighted Regression Adjustment." arXiv:2208.01300.

## Installation

```r
# install.packages("devtools")
devtools::install_github("deryauysal/drlate", subdir = "R")
```

## Example

```r
library(drlate)
sipp <- haven::read_dta(system.file("extdata", "sipp.dta", package = "drlate"))
?drlate   # see the examples in the help file
```

## License

MIT, see `LICENSE.md`.
