## A63: NHN_Z/NE_Z fits whose sigma_v reaches its lower bound used to report
## standard errors from a Hessian that has no meaning there.

test_that(".sigma_v_corner() applies NTN's 1e-3 rule relative to sigma_u", {
  expect_true(.sigma_v_corner(1e-7, 1))
  expect_true(.sigma_v_corner(9e-4, 1))
  expect_false(.sigma_v_corner(1.1e-3, 1))
  expect_false(.sigma_v_corner(0.3, 1))
  expect_false(.sigma_v_corner(NA_real_, 1))
})

test_that("a sigma_v boundary fit is flagged, warned, and carries no standard errors", {
  skip_on_cran()
  ## Whichever side of the bound each fit lands on, the report must match it.
  for (mn in c("NHN_Z", "NE_Z")) {
    for (s in c(2, 6, 42)) {
      d <- as.data.frame(data_gen_cs(N = 200, rand = s, sig_u = 1, sig_v = 0.3,
        cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5))
      w <- character()
      f <- withCallingHandlers(
        sfm(y_pcs_z ~ x1 + x2 | z, model_name = mn, data = d),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        }
      )
      info <- paste(mn, "seed", s)
      expect_type(f$sigma_v_at_bound, "logical")
      if (isTRUE(f$sigma_v_at_bound)) {
        expect_true(all(is.na(f$out[, "st_err"])), info = info)
        expect_true(any(grepl("sigma_v has collapsed to its lower bound", w)), info = info)
        expect_false(any(grepl("not a maximum", w)), info = info)
      } else {
        expect_false(any(grepl("sigma_v has collapsed", w)), info = info)
      }
    }
  }
})
