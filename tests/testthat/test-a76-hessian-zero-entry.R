## A76. psfm() discarded every standard error when ANY Hessian entry was
## exactly zero. A scale flat at its bound has exact-zero cross-derivatives,
## so on Ubuntu R oldrel the GTRE_FML fit on the A55 design lost all of its
## standard errors, including those with clearly positive variances.

test_that("one zero entry does not make a Hessian unusable", {
  h <- matrix(c(4, 0, 1, 0, 3, 0, 1, 0, 2), 3, 3)
  expect_false(sfa:::.hessian_unusable(h))
  expect_true(all(is.finite(sqrt(diag(solve(h))))))
})

test_that("the placeholders that do mean no Hessian are caught", {
  expect_true(sfa:::.hessian_unusable(NULL))
  expect_true(sfa:::.hessian_unusable(matrix(NA_real_, 2, 2)))
  expect_true(sfa:::.hessian_unusable(matrix(0, 3, 3)))
  expect_true(sfa:::.hessian_unusable(matrix(c(1, Inf, Inf, 1), 2, 2)))
  expect_true(sfa:::.hessian_unusable(1))
})

test_that("psfm() no longer tests for a zero entry", {
  src <- paste(deparse(body(psfm)), collapse = "\n")
  expect_false(grepl("any\\((opt|fit)\\$hessian == 0", src))
})
