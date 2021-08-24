# gKRLS

Chang and Goplerud (n.d.)'s implementation of KRLS as a mixed effect model. This package modifies the underlying architecture of the Laplace approximations in `lme4` (link; following Ben Bolker's [notes](https://bbolker.github.io/mixedmodels-misc/notes/varmats.html)) to fit kernel regressions for any non-linear outcome that can be used in `glmer`. We also provide an implementation of random sketching or projection following [Yang et al. (2017)](https://doi.org/10.1214/16-AOS1472).

The syntax is straightforward to users of KRLS. Specifically, the following command will fit our generalized KRLS algorithm. Fixed effects can be included alongside the kernel using the `formula` argument; variables that go into the kernel are put in `kernel_X`. The following example estimates a Poisson regression with an intercept and a flexible kernel term.

```
gKRLS(formula = y ~ 1,  kernel_X = X_train, data = data, family = poisson())
 ```

Sketching is automatically applied such that the dimensionality of the sketched problem is `5 * ceiling(N^(1/3))`. This can be modified directly by the user with the `sketch_dimension = X` arguments. Setting this to `NULL` results in no sketching. Eigentruncation can also be employed by setting the `eigentruncation = X` to some tolerance. This is not used by default but if used in conjunction with sketching, it operates on the eigenvalues of the sketched matrix.
