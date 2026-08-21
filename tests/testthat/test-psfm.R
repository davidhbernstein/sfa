## Panel entry point.

test_that("every psfm() model_name fits and returns a well-formed sfareg", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 50)
  specs <- list(
    list("TFE",           y_tfe  ~ x1_w + x2_w),
    list("TFE_WMLE",      y_tfe  ~ x1_w + x2_w),
    list("SSFE",          y_ssfe ~ x1   + x2),
    list("FD",            y_fd   ~ x_fd | z_fd),
    list("TRE",           y_tre  ~ x1   + x2),
    list("GTRE",          y_gtre ~ x1   + x2),
    list("PL80",          y_ssfe ~ x1   + x2),
    list("BC92",          y_bc92 ~ x1   + x2),
    list("K1990",         y_bc92 ~ x1   + x2),
    list("K1990modified", y_bc92 ~ x1   + x2)
  )
  for (s in specs) {
    fit <- fit_tfe_quietly(psfm(s[[2]], model_name = s[[1]], data = d, individual = "name"))
    expect_s3_class(fit, "sfareg")
    expect_identical(fit$model_name, s[[1]])
    expect_true(all(is.finite(fit$coefficients)), info = s[[1]])
    expect_equal(unname(fit$coefficients), unname(fit$out[, "par"]), info = s[[1]])
  }
})

test_that("psfm() accepts a plain data.frame as well as a pdata.frame", {
  skip_on_cran()
  skip_if_not_installed("plm")
  d  <- panel_small()
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  a  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d,  individual = "name")
  b  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = pd, individual = "name")
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = 1e-6)
})

test_that("PL80 and BC92 recover their true parameters", {
  skip_on_cran()
  ## y_ssfe uses time-invariant inefficiency (sig_u = 1), which is exactly
  ## PL80's assumption; y_bc92 adds eta = 0.1 decay.
  d <- panel_small(t = 8, N = 150)
  p <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d, individual = "name")
  expect_equal(unname(p$coefficients["x1"]), 0.5, tolerance = 0.05)
  expect_equal(unname(p$coefficients["x2"]), 0.5, tolerance = 0.05)

  b <- psfm(y_bc92 ~ x1 + x2, model_name = "BC92", data = d, individual = "name")
  expect_equal(unname(b$coefficients["x1"]), 0.5, tolerance = 0.05)
  ## The decay parameter is reported as "time" (frontier's naming); the
  ## generator uses eta = 0.1.
  expect_equal(unname(b$coefficients["time"]), 0.1, tolerance = 0.06)
})

test_that("PL80 is the eta = 0 special case of BC92", {
  skip_on_cran()
  ## Both share one likelihood, differing only in B_it. On a time-invariant
  ## DGP, BC92's eta should be near zero and its fit close to PL80's.
  d <- panel_small(t = 8, N = 120)
  p <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d, individual = "name")
  b <- psfm(y_ssfe ~ x1 + x2, model_name = "BC92", data = d, individual = "name")
  expect_lt(abs(unname(b$coefficients["time"])), 0.05)
  expect_equal(unname(b$coefficients["x1"]), unname(p$coefficients["x1"]), tolerance = 0.02)
  ## BC92 nests PL80, so it cannot fit worse.
  expect_gte(as.numeric(logLik(b)), as.numeric(logLik(p)) - 1e-4)
})

test_that("K1990 and K1990modified nest PL80's constant path", {
  skip_on_cran()
  d  <- panel_small(t = 8, N = 100)
  p  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80",  data = d, individual = "name")
  for (mn in c("K1990", "K1990modified")) {
    k <- psfm(y_ssfe ~ x1 + x2, model_name = mn, data = d, individual = "name")
    expect_true(all(is.finite(k$coefficients)), info = mn)
    expect_gte(as.numeric(logLik(k)), as.numeric(logLik(p)) - 1e-4)
  }
})

test_that("SSFE is deterministic and reports no log-likelihood", {
  d   <- panel_small()
  fit <- psfm(y_ssfe ~ x1 + x2, model_name = "SSFE", data = d, individual = "name")
  expect_null(fit$opt)
  ## Not an MLE, so logLik()/AIC()/BIC() must decline rather than invent a value.
  expect_warning(ll <- logLik(fit))
  expect_true(is.na(ll))
  expect_length(fit$alpha_hat, length(unique(d$name)))
  expect_true(all(fit$exp_u_hat > 0 & fit$exp_u_hat <= 1))
  ## Deterministic: refitting gives bit-identical numbers.
  again <- psfm(y_ssfe ~ x1 + x2, model_name = "SSFE", data = d, individual = "name")
  expect_identical(fit$coefficients, again$coefficients)
})

test_that("psfm() reports collinear between-individual designs rather than failing opaquely", {
  skip_on_cran()
  ## Time dummies are estimable in the pooled design but collapse onto the
  ## intercept once averaged within each individual, which is what breaks the
  ## random-effects starting-value regression.
  d          <- panel_small(t = 5, N = 40)
  d$year_fac <- factor(d$year)
  f          <- y_tre ~ x1 + x2 + year_fac

  expect_error(psfm(f, model_name = "TRE", data = d, individual = "name",
                    collinear_action = "error"), "collinear|rank")
  expect_warning(psfm(f, model_name = "TRE", data = d, individual = "name",
                      collinear_action = "start_only"))
})

test_that("psfm() rejects an unknown model_name", {
  d <- panel_small()
  expect_error(psfm(y_tre ~ x1 + x2, model_name = "NOPE", data = d, individual = "name"))
})

test_that("pipe-count validation rejects a formula the model cannot use", {
  d <- panel_small()
  ## TFE takes no pipe segments; TRE requires the explicit _Z name to use one.
  expect_error(psfm(y_tfe ~ x1_w + x2_w | z_fd, model_name = "TFE",
                    data = d, individual = "name"))
  expect_error(psfm(y_tre ~ x1 + x2 | z_gtre, model_name = "TRE",
                    data = d, individual = "name"))
})
