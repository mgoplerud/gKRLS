# gKRLS 
[![CRAN status](https://www.r-pkg.org/badges/version/gKRLS)](https://CRAN.R-project.org/package=gKRLS) [![R-CMD-check](https://github.com/mgoplerud/gKRLS/workflows/R-CMD-check/badge.svg)](https://github.com/mgoplerud/gKRLS/actions) [![codecov](https://codecov.io/gh/mgoplerud/gKRLS/branch/cran/graph/badge.svg?token=U22YCB3LPU)](https://app.codecov.io/gh/mgoplerud/gKRLS)

This package implements [Chang and Goplerud (2023)](https://doi.org/10.1017/pan.2023.27)'s generalization of Kernel Regularized Least Squares (gKRLS), also known as kernel ridge regression. This reformulates [g]KRLS as a hierarchical model. Estimation proceeds using `mgcv` and associated functions such as `gam`, `bam`, or `gamm4`. Thus, one can use `gKRLS` for any outcome implemented in `mgcv` as well as including multiple smooth terms, non-penalized covariates, etc. We also provide an implementation of random sketching following [Yang et al. (2017)](https://doi.org/10.1214/16-AOS1472).

The package can be installed from CRAN or the most-to-update version can be installed using `devtools`.

```
# CRAN
install.packages("gKRLS")
# Up-to-Date GitHub Version
library(remotes)
remotes::install_github("mgoplerud/gKRLS", dependencies = TRUE)
```

The syntax is straightforward to users of `mgcv`. The following example estimates a Poisson regression with an intercept and a flexible kernel term.

```
gam(y ~ s(x1, x2, bs = "gKRLS"), data = data, family = poisson())
 ```

`gKRLS` by default uses subsampling sketching (i.e., building the kernel based on a random sample of observations) where the dimensionality of the sketched kernel is `5 * ceiling(N^(1/3))`. Using `xt = gKRLS(...)` can modify the type of sketching. Please see the documentation for details.

Functions are also available to implement `gKRLS` in an ensemble using `SuperLearner` and in double/debiased machine learning using `DoubleML`. It also allows `sandwich` to calculate robust or clustered standard errors for standard families when using `gam` or `bam`; see [Chang and Goplerud (2023)](https://doi.org/10.1017/pan.2023.27) for more details.

`calculate_effects` can compute average marginal effects and predicted values. The examples for `calculate_effects` show how to calculate quantities such as predicted probability curves.
