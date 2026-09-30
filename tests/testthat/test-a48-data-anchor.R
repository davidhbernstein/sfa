## A48. .sfa_data_check() validates recovered data in tiers. sfm() stores its
## OLS residuals, so a recovered frame is rebuilt and compared -- decisive.
## zsfm() and ttsfm() stored no residuals: ZISF, TTNE and TTHN reached only
## the row-count tier, and TTNLS stored nothing at all. A decoy with the same
## rows and columns, bound to the name the call recorded, passed the check,
## and fitted()/residuals()/predict() answered from it without a word.

.a48_data <- function(seed) {
  data_gen_cs(N = 120, rand = seed, sig_u = 0.5, sig_v = 0.2, cons = 5,
    beta1 = 1.5, beta2 = 2.0, a = 5, mu = 0.1)
}

## Each call names `a48_d` itself, so the call the fit records is the one
## rebound below. (Fitting through a wrapper would record the wrapper's own
## argument instead, and the decoy would never be looked up.)
.a48_fits <- list(
  ZISF  = quote(zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = a48_d)),
  TTNE  = quote(ttsfm(y_ttne ~ x1 + x2, model_name = "TTNE", data = a48_d)),
  TTHN  = quote(ttsfm(y_tthn ~ x1 + x2, model_name = "TTHN", data = a48_d)),
  TTNLS = quote(ttsfm(y_ttne ~ x1 + x2, model_name = "TTNLS", data = a48_d))
)

## Fit on `a48_d` in a fresh environment, so the name can be rebound there.
.a48_fit_in <- function(model) {
  e <- new.env(parent = globalenv())
  e$a48_d <- .a48_data(42)
  e$fit <- suppressWarnings(eval(.a48_fits[[model]], e))
  e
}

for (m in names(.a48_fits)) {
  test_that(paste(m, "stores an anchor that matches its own data"), {
    skip_on_cran()
    e <- .a48_fit_in(m)
    expect_false(is.null(e$fit$anchor_resid))
    expect_identical(e$fit$nobs, nrow(e$a48_d))
    ## Positive control: the stored residuals are exactly what the check
    ## rebuilds from the real data, so the real data still passes.
    chk <- .sfa_data_check(e$fit, e$a48_d)
    expect_identical(chk$checked, "ols_residuals")
    expect_true(chk$ok)
    ## And the private name leaves the ols_residuals consumers alone.
    expect_null(e$fit$ols_residuals)
  })

  test_that(paste(m, "refuses a same-shaped decoy instead of answering from it"), {
    skip_on_cran()
    e <- .a48_fit_in(m)
    decoy <- .a48_data(999)
    expect_identical(dim(decoy), dim(e$a48_d))
    expect_false(isTRUE(all.equal(decoy$x1, e$a48_d$x1)))
    e$a48_d <- decoy
    expect_false(.sfa_data_check(e$fit, decoy)$ok)
    expect_error(stats::fitted(e$fit), "not the data this model was fitted to")
  })
}
