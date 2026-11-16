## Two reporting defects from the 2026-11-15 sweep over 42 model
## specifications. Both are asserted against objects BUILT BY HAND: what is
## under test is what the code says about a fit, not where an optimizer lands,
## and a hand-built object pins the arithmetic without needing either.

## ---- A81: plot.sfareg(which = 1) on an all-zero Hessian --------------------
## Panel 1 draws the eigenvalue spectrum with log = "y". A fit whose stored
## Hessian is all zeros has every eigenvalue 0, so there is no positive range,
## and plot() used to STOP -- leaving the diagnostic unavailable on exactly the
## fits whose diagnostics matter most.

zero_hess_fit <- function(p = 4L) {
  nm <- c("(Intercept)", paste0("x", seq_len(p - 1L)))
  out <- cbind(par = rep(0.5, p), st_err = rep(NA_real_, p), `t-val` = rep(NA_real_, p))
  rownames(out) <- nm
  z <- list(
    model_name = "NHN", out = out,
    coefficients = stats::setNames(out[, "par"], nm),
    std.errors = stats::setNames(out[, "st_err"], nm),
    opt = list(hessian = matrix(0, p, p, dimnames = list(nm, nm)), value = 10),
    nobs = 100L
  )
  class(z) <- "sfareg"
  z
}

test_that("the Hessian-spectrum panel draws on an all-zero Hessian", {
  f <- zero_hess_fit()
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  ## The assertion is that it DRAWS. Before the fix this raised from plot().
  expect_error(plot(f, which = 1L), NA)
})

test_that("a spectrum with some positive eigenvalues still uses the log panel", {
  f <- zero_hess_fit()
  h <- f$opt$hessian
  diag(h) <- c(1, 2, 0, 0)          # two positive, two not
  f$opt$hessian <- h
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_error(plot(f, which = 1L), NA)
})

## ---- A82: print.sfma()'s model count ---------------------------------------
## The line was built from length(x$models) on BOTH sides of the "of", so a run
## that dropped a specification printed as though nothing had been dropped.

fake_sfma <- function(fitted, requested) {
  z <- list(
    weights = stats::setNames(rep(1 / length(fitted), length(fitted)), fitted),
    u_hat = rep(0.3, 5), efficiency = rep(exp(-0.3), 5),
    models = fitted, models_requested = requested,
    k = rep(4L, length(fitted)),
    logLik = rep(-100, length(fitted)), nobs = 5L,
    weight_type = "sfma",
    ## print.sfma() reads criterion$full_model and u_by_model, so the fixture
    ## has to carry the shape the real object has.
    criterion = list(full_model = fitted[1], value = 200),
    u_by_model = matrix(0.3, nrow = 5, ncol = length(fitted),
                        dimnames = list(NULL, fitted))
  )
  class(z) <- "sfma"
  z
}

test_that("print.sfma reports the requested count, not the surviving one twice", {
  z <- fake_sfma(c("NHN", "NE"), c("NHN", "NE", "NTN"))
  out <- capture.output(print(z))
  expect_true(any(grepl("2 of 3 fitted", out, fixed = TRUE)))
  expect_false(any(grepl("3 of 3 fitted", out, fixed = TRUE)))
  ## and it names what went missing, which is the point of the line
  expect_true(any(grepl("dropped", out) & grepl("NTN", out)))
})

test_that("print.sfma says nothing about drops when nothing was dropped", {
  z <- fake_sfma(c("NHN", "NE", "NTN"), c("NHN", "NE", "NTN"))
  out <- capture.output(print(z))
  expect_true(any(grepl("3 of 3 fitted", out, fixed = TRUE)))
  expect_false(any(grepl("dropped", out)))
})

test_that("an sfma object from before the fix does not get an invented denominator", {
  z <- fake_sfma(c("NHN", "NE"), c("NHN", "NE", "NTN"))
  z$models_requested <- NULL
  out <- capture.output(print(z))
  expect_true(any(grepl("2 fitted", out, fixed = TRUE)))
  expect_false(any(grepl(" of ", out, fixed = TRUE)))
})
