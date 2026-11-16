## A77: sfm() nulls every standard error at the deterministic-frontier corner,
## but vcov() went back to the stored Hessian -- which still inverts there --
## and handed the withheld numbers out again. lmtest::coeftest() takes its
## standard errors from vcov(), not from $out, so it reported z = 1.97e6 on a
## fit whose own summary() printed NA for the same coefficient.
##
## The fixtures below are built by hand rather than fitted: the contract is
## "whatever the fit reports as NA, vcov() reports as NA", and asserting it
## must not depend on an optimizer reaching the corner.

mk_fit <- function(se, flag = NULL) {
  nm <- c("sigma_v", "(Intercept)", "x1")
  p <- length(nm)
  ## A well-conditioned Hessian, so the default path definitely inverts and
  ## the mask is the only thing that can produce an NA.
  H <- diag(c(4, 2, 8))
  H[1, 2] <- H[2, 1] <- 0.5
  out <- cbind(par = c(0.3, 0.5, 0.5), st_err = se, `t-val` = c(1, 1, 1) / se)
  rownames(out) <- nm
  f <- list(
    coefficients = stats::setNames(c(0.3, 0.5, 0.5), nm),
    std.errors = stats::setNames(se, nm),
    out = out,
    opt = list(hessian = H, par = c(0.3, 0.5, 0.5)),
    par_index = NULL, par_scale = NULL
  )
  if (!is.null(flag)) f[[flag]] <- TRUE
  class(f) <- "sfareg"
  f
}

test_that("a boundary-corner fit gets NAs from vcov(), not the inverted Hessian", {
  for (flag in c("sigma_v_at_bound", "ntn_sigma_v_at_bound")) {
    f <- mk_fit(rep(NA_real_, 3), flag)
    ## The Hessian on this object inverts cleanly, which is the whole point:
    ## invertibility is not evidence that a standard error exists.
    expect_true(all(is.finite(solve(f$opt$hessian))))
    expect_warning(V <- vcov(f), "no finite standard error exists")
    expect_true(all(is.na(V)), info = flag)
    expect_identical(dim(V), c(3L, 3L))
    expect_identical(rownames(V), names(coef(f)))
    ## The scope, pinned: the check belongs to the DEFAULT path, which is the
    ## one coeftest(), confint() and sandwich::sandwich() reach. A type asked
    ## for by name is a different estimator and is left to its own behaviour
    ## -- here, its own error about the missing score matrix, which proves the
    ## boundary check did not intercept it. test-vcov-bhhh-advice.R owns the
    ## contract that those types error informatively on such a fit.
    expect_error(vcov(f, type = "bhhh"), "score matrix is unavailable")
  }
})

test_that("vcov() reports NA exactly where the fit reports NA", {
  ## One unresolved standard error, not a whole corner: A75's case.
  f <- mk_fit(c(0.1, NA_real_, 0.2))
  expect_warning(V <- vcov(f), "reports no standard error for \\(Intercept\\)")
  expect_true(all(is.na(V[2, ])))
  expect_true(all(is.na(V[, 2])))
  expect_true(all(is.finite(V[c(1, 3), c(1, 3)])))
  ## The invariant, stated positionally -- the two identical "(Intercept)"
  ## rows A80 leaves on some fits make a name-based comparison self-match.
  se_fit <- f$out[, "st_err"]
  se_vcov <- suppressWarnings(sqrt(diag(V)))
  expect_identical(which(is.na(se_fit)), which(is.na(se_vcov)))
  ## And the standard errors that DO survive still agree with the Hessian.
  expect_equal(unname(V[1, 1]), solve(f$opt$hessian)[1, 1], tolerance = 1e-12)
})

test_that("a fit with every standard error reported is untouched and silent", {
  f <- mk_fit(c(0.1, 0.2, 0.3))
  expect_silent(V <- vcov(f))
  expect_equal(unname(V), solve(f$opt$hessian), tolerance = 1e-12)
})

test_that("the diagonal fallback masks a withheld parameter too", {
  ## No Hessian, so vcov() builds a diagonal from the reported standard
  ## errors. The withheld parameter's zero off-diagonals are as much of a
  ## claim as a non-zero one would be.
  f <- mk_fit(c(0.1, NA_real_, 0.2))
  f$opt$hessian <- NULL
  expect_warning(V <- vcov(f), "reports no standard error for")
  expect_true(all(is.na(V[2, ])))
  expect_equal(unname(V[1, 1]), 0.01, tolerance = 1e-12)
})

test_that("on a real sigma_v-corner fit, vcov() and $out agree", {
  skip_on_cran()
  ## Whichever side of the bound each fit lands on, vcov() must match the
  ## report -- so this asserts the agreement, never that the corner is reached.
  for (mn in c("NHN_Z", "NE_Z")) {
    for (s in c(2, 7, 42)) {
      d <- as.data.frame(data_gen_cs(N = 200, rand = s, sig_u = 1, sig_v = 0.3,
        cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5))
      f <- suppressWarnings(
        sfm(y_pcs_z ~ x1 + x2 | z, model_name = mn, data = d)
      )
      info <- paste(mn, "seed", s)
      V <- suppressWarnings(vcov(f))
      expect_identical(
        which(is.na(f$out[, "st_err"])),
        which(is.na(suppressWarnings(sqrt(diag(V))))),
        info = info
      )
      if (isTRUE(f$sigma_v_at_bound)) {
        expect_true(all(is.na(V)), info = info)
        ## The accessor the defect was measured through.
        skip_if_not_installed("lmtest")
        ct <- suppressWarnings(lmtest::coeftest(f))
        expect_true(all(is.na(ct[, 2])), info = info)
      }
    }
  }
})
