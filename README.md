# gKRLS [![R-CMD-check](https://github.com/mgoplerud/gKRLS/workflows/R-CMD-check/badge.svg)](https://github.com/mgoplerud/gKRLS/actions) [![Codecov test coverage](https://codecov.io/gh/mgoplerud/gKRLS/branch/cran/graph/badge.svg)](https://app.codecov.io/gh/mgoplerud/gKRLS?branch=cran)

This package implements Chang and Goplerud (2022)'s generalization of Kernel Regularized Least Squares (gKRLS), also known as kernel ridge regression. This reformulates [g]KRLS as a hierarchical model. Estimation proceeds using \code{mgcv} and associated functions such as \code{gam}, \code{bam}, or \code{gamm4}. Thus, it can be used for any outcome implemented in that software as well as including multiple smooth terms, non-penalized covariates, etc.

We also provide an implementation of random sketching or projection following [Yang et al. (2017)](https://doi.org/10.1214/16-AOS1472).

The syntax is straightforward to users of \texttt{mgcv}. The following example estimates a Poisson regression with an intercept and a flexible kernel term.

```
gam(y ~ s(x1, x2, bs = "gKRLS"), data = data, family = poisson())
 ```

Sketching is automatically applied such that the dimensionality of the sketched problem is `5 * ceiling(N^(1/3))`. This can be modified directly by the user with the `xt = gKRLS(...)` arguments.

Functions are also available to implement `gKRLS` in an ensemble using `SuperLearner` and in debiased/double machine learning using `DoubleML'.

Marginal effects can be calculated using the `calculate_effects` function. Please see the documentation for details.