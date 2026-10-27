## A65. zsfm(), lcsfm(), copsfm(), ivsfm() and selsfm() reported a NaN standard
## error with no signal whatsoever: suppressWarnings() discarded R's own "NaNs
## produced", which was the only thing said about it, nothing set a flag, and
## the convergence code read 0. They now use .se_missing_report(), the
## model-agnostic helper PR #63 added for ttsfm(), and guard solve() so a
## singular-but-nonzero Hessian stops discarding a finished fit.
##
## Whether a given design LANDS on an indefinite Hessian is platform-dependent
## -- the optimizer's stopping point, and the sign of a small eigenvalue there,
## is BLAS-dependent arithmetic -- so nothing below pins which parameters lose
## a standard error, or how many. Every expectation is derived from the fit in
## front of it, and the regime itself is a skip condition rather than an
## assertion. The helper's own contract is unit-tested deterministically in
## test-a53-missing-se-report.R, where the input is constructed.

.a65_zisf <- function(N, seed) {
  as.data.frame(data_gen_cs(N = N, rand = seed, sig_u = 1, sig_v = 0.3,
                            cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                            mu = 0.5, a = 5))
}

test_that("zsfm('ZISF') names the standard errors it could not compute", {
  skip_on_cran()
  ## N = 25, seed 8. On macOS arm64, where this was found, the Hessian at the
  ## reported point has one negative eigenvalue (-1474 against a largest of
  ## 3.48e6), two of six standard errors go missing -- one of them the slope
  ## x2, whose true value is 0.5 -- and conv is 0, so nothing else in the
  ## output flags it. Half a second to fit.
  w <- character(0)
  f <- withCallingHandlers(
    zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = .a65_zisf(25, 8)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  se <- f$out[, "st_err"]
  bad <- names(se)[!is.finite(se)]
  skip_if_not(length(bad) > 0L,
    "every standard error is finite on this platform, so there is nothing to report"
  )

  r <- f$se_unresolved
  expect_false(is.null(r))
  ## Named for exactly the parameters that lost one, however many that is here.
  expect_setequal(r$parameters, bad)
  expect_identical(r$n_missing, length(bad))
  expect_identical(r$n_par, length(se))
  ## The warning fired, and it is the one this change adds.
  expect_true(any(grepl("no standard error could be computed", w, fixed = TRUE)))
  expect_match(r$message, "zsfm(model_name = \"ZISF\")", fixed = TRUE)
  ## The two readings are not distinguished, and the message says so.
  expect_match(r$message, "EITHER a point that is not a maximum OR a boundary")
  ## Internally consistent: "indefinite" iff a negative eigenvalue was read.
  expect_identical(isTRUE(r$indefinite), isTRUE(r$min_eigenvalue < 0))

  ## The flag, not just the warning: a warning does not survive saveRDS().
  out <- paste(capture.output(print(f)), collapse = " ")
  expect_match(out, "no standard error could be computed for", fixed = TRUE)
  expect_match(out, paste0(length(bad), " of ", length(se), " parameters"),
    fixed = TRUE
  )
})

test_that("an ordinary ZISF fit is left alone", {
  skip_on_cran()
  ## The same design at N = 300: on this Mac the min eigenvalue is 37.1 and all
  ## six standard errors are finite. The claim being tested is that a HEALTHY
  ## fit gets no flag, so a platform whose fit is not healthy skips instead of
  ## failing -- that case is the test above, not this one.
  f <- suppressWarnings(
    zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = .a65_zisf(300, 1))
  )
  skip_if_not(all(is.finite(f$out[, "st_err"])),
    "this platform's N = 300 fit lost a standard error of its own"
  )
  expect_null(f$se_unresolved)
  out <- paste(capture.output(print(f)), collapse = " ")
  expect_false(grepl("no standard error could be computed", out, fixed = TRUE))
})
