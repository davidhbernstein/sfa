## A persistent scale on the zero boundary in a GTRE fit.

test_that("the helper judges a collapse relative to scale", {
  expect_true(sfa:::.panel_scale_at_bound(1e-7, 1.0))
  expect_false(sfa:::.panel_scale_at_bound(0.4, 1.0))
  ## Relative, not absolute: 0.005 is a collapse against a scale of 1 but not
  ## against a scale of 0.05.
  expect_true(sfa:::.panel_scale_at_bound(0.005, 1.0))
  expect_false(sfa:::.panel_scale_at_bound(0.005, 0.05))
  ## Degenerate inputs return NA rather than a confident answer.
  expect_true(is.na(sfa:::.panel_scale_at_bound(NA_real_, 1)))
  expect_true(is.na(sfa:::.panel_scale_at_bound(0.1, 0)))
  expect_true(is.na(sfa:::.panel_scale_at_bound(0.1, NA_real_)))
})

test_that("GTRE reports a collapsed persistent scale, in both directions", {
  skip_on_cran()
  ## Two seeds that collapse in OPPOSITE directions -- 1005 sends sigh to the
  ## boundary, 1001 sends sigr. Both are genuine maximum likelihood solutions:
  ## on seed 1005 the boundary fit beats the TRUE parameter vector by 3.66 log
  ## units. See notes/sigh_investigation_2026-08-29.md.
  grab <- function(seed) {
    d <- data_gen_p(t = 6, N = 200, rand = seed, sig_u = 1, sig_v = 0.3,
                    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5,
                    beta2 = 0.5)
    pd <- plm::pdata.frame(d, index = c("name", "year"))
    w <- character(0)
    f <- withCallingHandlers(
      psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name",
           halton_num = 50, rand.gtre = 7L, estimator = "sml"),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    list(fit = f, warns = w)
  }

  a <- grab(1005)
  expect_true(isTRUE(a$fit$sigh_at_bound))
  expect_false(isTRUE(a$fit$sigr_at_bound))
  expect_true(any(grepl("sigh has collapsed", a$warns)))
  ## The warning must name the OTHER scale, because they have to be read
  ## together -- the variation is reclassified, not lost.
  expect_true(any(grepl("sigr = ", a$warns)))

  b <- grab(1001)
  expect_true(isTRUE(b$fit$sigr_at_bound))
  expect_false(isTRUE(b$fit$sigh_at_bound))
  expect_true(any(grepl("sigr has collapsed", b$warns)))
  expect_true(any(grepl("sigh = ", b$warns)))
})

test_that("the warning explains that the boundary can be the correct MLE", {
  skip_on_cran()
  d <- data_gen_p(t = 6, N = 200, rand = 1005, sig_u = 1, sig_v = 0.3,
                  sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  expect_warning(
    psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name", halton_num = 50,
         rand.gtre = 7L, estimator = "sml"),
    "CORRECT maximum likelihood estimate"
  )
})

test_that("psfm_bootstrap warns when bootstrapping a boundary fit", {
  skip_on_cran()
  ## A parametric bootstrap resamples from the FITTED model, so a collapsed
  ## persistent scale means resampling from a DGP without that component --
  ## and the bootstrap is inconsistent on the boundary of the parameter space
  ## regardless. The interval would understate the uncertainty in the
  ## persistent split rather than represent it.
  d <- data_gen_p(t = 6, N = 200, rand = 1005, sig_u = 1, sig_v = 0.3,
                  sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  f <- suppressWarnings(
    psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name", halton_num = 50,
         rand.gtre = 7L, estimator = "sml"))
  expect_true(isTRUE(f$sigh_at_bound))

  w <- character(0)
  invisible(withCallingHandlers(
    try(psfm_bootstrap(f, B = 2, numCores = 1, individual = "name"),
        silent = TRUE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }))
  expect_true(any(grepl("persistent scale on the zero boundary", w)))
})
