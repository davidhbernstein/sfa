## extract() for "sfareg", so texreg renders these fits into paper tables.

skip_if_no_texreg <- function() {
  testthat::skip_if_not_installed("texreg")
}

tex_fit <- function(seed = 1, n = 300) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.5 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n))
  sfm(y ~ x1 + x2, model_name = "NHN", data = data.frame(y, x1, x2))
}

test_that("extract() produces a texreg object carrying the fit's numbers", {
  skip_if_no_texreg()
  skip_on_cran()
  f <- tex_fit()
  tr <- extract.sfareg(f)

  expect_s4_class(tr, "texreg")
  expect_equal(tr@coef, unname(f$out[, "par"]))
  expect_equal(tr@se, unname(f$out[, "st_err"]))
  expect_identical(tr@coef.names, rownames(f$out))

  ## p-values are recomputed from par/se rather than read off the t-value
  ## column, so a fit with NA standard errors gives NA p-values instead of
  ## whatever happened to be stored.
  expect_equal(
    tr@pvalues,
    unname(2 * pnorm(-abs(f$out[, "par"] / f$out[, "st_err"])))
  )
  expect_true(all(c("Log Likelihood", "AIC", "BIC", "Num. obs.") %in% tr@gof.names))
  expect_equal(tr@gof[tr@gof.names == "Num. obs."], nobs(f))
})

test_that("include.scales separates the frontier from the variance parameters", {
  skip_if_no_texreg()
  skip_on_cran()
  f <- tex_fit()
  full <- extract.sfareg(f)
  slim <- extract.sfareg(f, include.scales = FALSE)

  expect_lt(length(slim@coef), length(full@coef))
  expect_true(all(c("x1", "x2") %in% slim@coef.names))
  ## Whatever the model's scale parameters are called, none should survive.
  expect_false(any(grepl("^(sigma|sig|lambda|gamma)", slim@coef.names)))
})

test_that("a moment-based fit renders without a log-likelihood row", {
  skip_if_no_texreg()
  skip_on_cran()
  ## C2SLS maximises nothing, so logLik() is NA. The row should be OMITTED
  ## rather than printed as NA, which would invite a comparison that cannot be
  ## made against a likelihood-based column beside it.
  set.seed(3)
  n <- 600
  w1 <- rnorm(n); w2 <- rnorm(n); x1 <- rnorm(n); eta <- rnorm(n)
  v <- 0.5 * (0.6 * eta + sqrt(1 - 0.36) * rnorm(n))
  x2 <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
  y <- 0.5 + 0.8 * x1 - 0.6 * x2 + v - abs(rnorm(n))
  d <- data.frame(y = y, x1 = x1, x2 = x2, w1 = w1, w2 = w2)

  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "C2SLS"
  )
  tr <- extract.sfareg(f)
  expect_false("Log Likelihood" %in% tr@gof.names)
  expect_false("AIC" %in% tr@gof.names)
  expect_true("Num. obs." %in% tr@gof.names)
  expect_silent(invisible(texreg::screenreg(list(tr))))
})

test_that("texreg DISPATCHES on a fit, not just on an extracted object", {
  ## This test used to grep NAMESPACE for the literal
  ## S3method(texreg::extract, sfareg) and assert the line was present. The
  ## line WAS present and registered nothing: texreg::extract is an S4
  ## generic, S3method() writes into the S3 methods table, and S4 dispatch
  ## never reads it (gap A79). So the assertion passed for a full release
  ## while texreg::screenreg(list(fit)) reported "Neither texreg nor broom
  ## supports models of class sfareg".
  ##
  ## The only form of this test that can fail when dispatch is broken is to
  ## hand screenreg() a FIT and look for the coefficients in the output.
  ## Everything else -- the NAMESPACE line, existsMethod(), or passing an
  ## already-extracted object as the tests above do -- passes either way.
  ## existsMethod("extract", "sfareg", where = asNamespace("texreg")) is in
  ## particular NOT a usable check: it reads FALSE on a build where dispatch
  ## demonstrably works, because .onLoad() registers against the generic
  ## rather than into texreg's own table.
  skip_if_no_texreg()
  f <- tex_fit()
  out <- capture.output(texreg::screenreg(list(f)))
  expect_true(any(grepl("x1", out, fixed = TRUE)))
  expect_true(any(grepl("x2", out, fixed = TRUE)))
  expect_false(any(grepl("does not support|Neither texreg", out)))
})
