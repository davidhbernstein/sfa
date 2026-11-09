## A75. sfm() could lose a standard error with no warning and no
## $se_unresolved: A65 wired the report into the other entry points, not this
## one. Which fits land on an indefinite Hessian is platform arithmetic, so the
## wiring is tested by forcing a missing SE, and the reproducers only check
## the invariant that a missing SE without a boundary flag is always reported.

test_that("sfm() reports a standard error it could not compute", {
  d <- as.data.frame(data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = .3, mu = .5, a = 5,
    cons = .5, beta1 = .5, beta2 = .5))
  local_mocked_bindings(
    .sfa_se_from_var = function(v) {
      v[!is.finite(v) | v < 0] <- NA_real_
      v[3] <- NA_real_
      sqrt(v)
    }
  )
  expect_warning(
    f <- sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN"),
    "sfm(model_name = \"NHN\"): no standard error could be computed for",
    fixed = TRUE
  )
  se <- f$out[, "st_err"]
  expect_identical(f$se_unresolved$parameters, names(se)[3])
  expect_identical(f$se_unresolved$n_par, length(se))
  out <- paste(capture.output(print(f)), collapse = " ")
  expect_match(out, "no standard error could be computed for", fixed = TRUE)
})

test_that("an ordinary sfm() fit carries no report", {
  d <- as.data.frame(data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = .3, mu = .5, a = 5,
    cons = .5, beta1 = .5, beta2 = .5))
  f <- sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN")
  expect_true(all(is.finite(f$out[, "st_err"])))
  expect_null(f$se_unresolved)
})

test_that("a missing NTN standard error without a boundary flag is reported", {
  skip_on_cran()
  ## The gap entry's reproducers: on macOS mu drifts toward its exponential
  ## limit and loses its SE. Whatever happens here, the invariant must hold.
  for (r in c(9L, 10L)) {
    d <- data_gen_cs(N = 150, rand = r, sig_u = 1, sig_v = .8, mu = .5,
      a = 5, cons = .5, beta1 = .5, beta2 = .5)
    w <- character(0)
    f <- withCallingHandlers(
      sfm(y_pcs_tn ~ x1 + x2, data = d, model_name = "NTN"),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    se <- f$out[, "st_err"]
    bad <- names(se)[!is.finite(se)]
    info <- paste("rand", r)
    if (length(bad) && !isTRUE(f$ntn_sigma_v_at_bound)) {
      expect_setequal(f$se_unresolved$parameters, bad)
      expect_true(any(grepl("no standard error could be computed", w,
        fixed = TRUE)), info = info)
    } else {
      expect_null(f$se_unresolved, info = info)
    }
  }
})
