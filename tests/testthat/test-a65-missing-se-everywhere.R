## A65. zsfm(), lcsfm(), copsfm(), ivsfm() and selsfm() reported a NaN standard
## error with no signal whatsoever: suppressWarnings() discarded R's own "NaNs
## produced", which was the only thing said about it, nothing set a flag, and
## the convergence code read 0. They now use .se_missing_report(), the
## model-agnostic helper PR #63 added for ttsfm(), and guard solve() so a
## singular-but-nonzero Hessian stops discarding a finished fit.

test_that("zsfm('ZISF') names the standard errors it could not compute", {
  skip_on_cran()
  ## N = 25, seed 8 of the ZISF design: the Hessian at the reported point has
  ## one negative eigenvalue (-1474 against a largest of 3.48e6), two of its
  ## six inverse-diagonal entries are negative, and the two standard errors
  ## that go missing include the slope x2, whose true value is 0.5. conv is 0,
  ## so nothing else in the output flags it. Fits in about half a second.
  d <- as.data.frame(data_gen_cs(N = 25, rand = 8, sig_u = 1, sig_v = 0.3,
                                 cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                                 mu = 0.5, a = 5))
  expect_warning(
    f <- zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d),
    "no standard error could be computed"
  )
  expect_identical(f$opt$convergence, 0L)

  r <- f$se_unresolved
  expect_false(is.null(r))
  expect_setequal(r$parameters, c("(Intercept)", "x2"))
  expect_identical(r$n_missing, 2L)
  expect_identical(r$n_par, 6L)
  expect_true(isTRUE(r$indefinite))
  expect_lt(r$min_eigenvalue, 0)
  expect_gt(r$max_eigenvalue, 0)
  ## The two readings are not distinguished, and the message says so.
  expect_match(r$message, "EITHER a point that is not a maximum OR a boundary")
  expect_match(r$message, "zsfm(model_name = \"ZISF\")", fixed = TRUE)

  ## The flag, not just the warning: a warning does not survive saveRDS().
  out <- paste(capture.output(print(f)), collapse = " ")
  expect_match(out, "no standard error could be computed for", fixed = TRUE)
  expect_match(out, "2 of 6 parameters", fixed = TRUE)
})

test_that("an ordinary ZISF fit is left alone", {
  skip_on_cran()
  ## The same design at N = 300: min eigenvalue 37.1, all six standard errors
  ## finite. No warning, no flag -- the report must not fire on a healthy fit.
  d <- as.data.frame(data_gen_cs(N = 300, rand = 1, sig_u = 1, sig_v = 0.3,
                                 cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                                 mu = 0.5, a = 5))
  expect_silent(f <- zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d))
  expect_null(f$se_unresolved)
  expect_true(all(is.finite(f$out[, "st_err"])))
  out <- paste(capture.output(print(f)), collapse = " ")
  expect_false(grepl("no standard error could be computed", out, fixed = TRUE))
})
