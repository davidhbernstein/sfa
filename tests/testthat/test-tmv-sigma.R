## Gap A38. tmvtnorm::ptmvnorm(), then used, stopped GTRE_Z's efficiency scores on the
## Windows and Ubuntu oldrel-1 runners with "sigma must be positive definite":
## a firm's posterior covariance had lost positive definiteness. .tmv_sigma()
## repairs a rounding-sized defect and refuses anything larger.

## tmvtnorm's own conditions (it is no longer a dependency, gap A25).
.tmv_check <- function(M) {
  if (!isSymmetric(M, tol = sqrt(.Machine$double.eps))) stop("sigma must be a symmetric matrix")
  if (any(diag(M) <= 0)) stop("sigma all diagonal elements must be positive")
  if (det(M) <= 0) stop("sigma must be positive definite")
  invisible(NULL)
}

.with_eigenvalues <- function(ev, seed) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(length(ev)^2), length(ev))))
  S <- Q %*% diag(ev) %*% t(Q)
  (S + t(S)) / 2
}

test_that(".tmv_sigma() returns a positive definite matrix unchanged", {
  S <- .with_eigenvalues(c(2.4, 1, 0.5, 0.1, 0.01, 0.0097), seed = 1)
  expect_null(.tmv_check(S))
  expect_identical(.tmv_sigma(S), S)
})

test_that(".tmv_sigma() repairs a rounding-sized loss of positive definiteness", {
  S <- .with_eigenvalues(c(2.4, 1, 0.5, 0.1, 0.01, -1e-9), seed = 2)
  expect_error(.tmv_check(S), "positive definite")
  R <- .tmv_sigma(S)
  expect_null(.tmv_check(R))
  expect_lt(max(abs(R - S)), 1e-8)
})

test_that(".tmv_sigma() refuses a covariance that is indefinite beyond rounding", {
  S <- .with_eigenvalues(c(2, 1, 0.5, 0.1, -0.3), seed = 3)
  expect_error(.tmv_sigma(S), "not positive definite")
  expect_error(.tmv_sigma(matrix(NA_real_, 3, 3)), "not finite")
})
