
X <- apply(X, MARGIN = 2, scale)
krlsout <- suppressWarnings(krls(X=X,y=y, print.level = 0))

krls_kernel <- krlsout$K
man_kernel <- create_sketched_kernel(X_test = X, X_train = X, tS = diag(nrow(X)), bandwidth = ncol(X))

stopifnot(max(abs(krls_kernel - man_kernel)) < 1e-7)


krls_KS <- krlsout$K %*% S
KS <- create_sketched_kernel(X_test = X, X_train = X, 
                             tS = t(S), bandwidth = ncol(X))

stopifnot(max(abs(krls_KS - KS)) < 1e-7)
