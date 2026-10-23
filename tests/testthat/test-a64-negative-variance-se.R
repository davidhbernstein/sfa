## A64. A negative variance estimate was clamped by pmax(v, 0) and reported as
## a standard error of exactly 0, with t = par/0 = Inf -- an infinitely precise
## estimate for a direction whose estimated variance was negative. The helper
## reports such a parameter as NA instead, and the end-to-end case is the fit
## the defect was found on.

test_that("a negative variance estimate is NA, not zero", {
  v <- c(a = 4, b = -0.000387, c = 0.25)
  se <- sfa:::.sfa_se_from_var(v)
  expect_equal(se[["a"]], 2)
  expect_equal(se[["c"]], 0.5)
  expect_true(is.na(se[["b"]]))
  ## The old expression returned 0 here, and par/0 is Inf.
  expect_false(isTRUE(se[["b"]] == 0))
  expect_identical(names(se), names(v))
})

test_that("non-finite variance entries are NA, and the rest survive", {
  se <- sfa:::.sfa_se_from_var(c(1, NA_real_, NaN, Inf, -Inf, 9))
  expect_equal(se[c(1, 6)], c(1, 3))
  expect_true(all(is.na(se[2:5])))
})

test_that("a zero variance is left at zero, which is a separate question", {
  ## Only the manufactured zeros are the defect; an exact zero that the
  ## algebra really produced is reported as it stands.
  expect_equal(sfa:::.sfa_se_from_var(0), 0)
})

test_that("a t-value built on the result is NA rather than Inf", {
  par <- c(x1 = 0.542644)
  se <- sfa:::.sfa_se_from_var(c(x1 = -3e-06))
  expect_true(is.na(par / se))
  expect_false(is.infinite(par / se))
})

test_that("psfm('GTRE_FML') reports no zero standard error on the A55 design", {
  skip_on_cran()
  ## Seed 1, N = 40, T = 5: three parameters (x1, sigv, sigh) have a NEGATIVE
  ## variance estimate here, and on main all three came back se = 0, t = Inf.
  d <- as.data.frame(data_gen_p(t = 5, N = 40, rand = 1, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5, eta = 0.1))
  f <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE_FML",
                             data = d, individual = "name"))
  se <- f$out[, "st_err"]
  tv <- f$out[, "t-val"]
  ## The variance estimates really are negative -- the fit is not being asked
  ## to produce a different Hessian, only to report it honestly.
  v <- diag(solve(f$opt$hessian))
  expect_gt(sum(v < 0), 0)
  expect_equal(sum(is.finite(se) & se == 0), 0L)
  expect_equal(sum(is.infinite(tv)), 0L)
  ## Exactly the negative-variance coordinates lose their standard error.
  expect_true(all(is.na(se[v < 0])))
  expect_true(all(is.finite(se[v > 0])))
})
