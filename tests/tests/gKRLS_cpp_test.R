
X <- matrix(rnorm(1000), ncol = 2)
X <- apply(X, MARGIN = 2, scale)
y <- rnorm(nrow(X))
krlsout <- suppressWarnings(KRLS::krls(X=X,y=y, print.level = 0))

krls_kernel <- krlsout$K
man_kernel <- gKRLS:::create_sketched_kernel(X_test = X, X_train = X, tS = diag(nrow(X)), bandwidth = ncol(X))

testthat::expect_equivalent(krls_kernel, man_kernel)


S <- matrix(rnorm(nrow(X) * 5), ncol = 5)
krls_KS <- krlsout$K %*% S
KS <- gKRLS:::create_sketched_kernel(X_test = X, X_train = X, 
                             tS = t(S), bandwidth = ncol(X))

testthat::expect_equivalent(krls_KS, KS)

stop()

m <- build_kern_parallel(X_train = X, X_test = X, tS = t(S), 
                                 bandwidth = ncol(X), grain_size = 100,
                                 threads = 2)
stopifnot(max(abs(KS - m)) < 1e-7)


stop()

out <- data.frame()
for (it in 1:1){
  print(it)
  X <- matrix(rnorm(10^5 * 2), ncol = 2)
  X <- apply(X, MARGIN = 2, scale)
  size <- ceiling(nrow(X)^(1/3)) * 5
  S <- matrix(rnorm(nrow(X) * size), ncol = size)
  print('parl')
  
  time_parallel3 <- proc.time()
  parallel_KS3 <- build_kern_parallel(X_train = X, X_test = X, tS = t(S), 
                                      bandwidth = ncol(X), grain_size = 100,
                                      threads = 3)
  time_parallel3 <- proc.time() - time_parallel3
  
  time_parallel <- proc.time()
  parallel_KS <- build_kern_parallel(X_train = X, X_test = X, tS = t(S), 
                                             bandwidth = ncol(X), grain_size = 100,
                                             threads = 1)
  time_parallel <- proc.time() - time_parallel

  gc()
  
  
  gc()
  
  print('non')
  time_KS <- proc.time()
  KS <- create_sketched_kernel(X_test = X, X_train = X, 
                               tS = t(S), bandwidth = ncol(X))
  time_KS <- proc.time() - time_KS
  
  gc()
  stopifnot(max(abs(KS - parallel_KS)) < 1e-7)
  stopifnot(max(abs(KS - parallel_KS3)) < 1e-7)
  out <- rbind(out, data.frame(parallel = time_parallel[3], 
                               parallel_3 = time_parallel3[3],
             ks = time_KS[3], seed = it))
  print(out)
}
