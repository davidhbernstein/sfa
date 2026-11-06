## A transformed response (log(y) ~ ...) must be fitted as transformed. Until
## 1.2.1 every entry point built on data_proc() took the response as
## all.vars(formula)[1] -- the raw column -- while the regressors went through
## model.frame(), so log(q) ~ log(K) fitted q in levels on log(K).

expect_same_fit <- function(a, b, tol = 1e-6) {
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = tol)
  expect_equal(-a$opt$value, -b$opt$value, tolerance = tol)
}

test_that(".response_name() keeps a transformed left-hand side as written", {
  expect_identical(.response_name(y ~ x), "y")
  expect_identical(.response_name(log(q) ~ log(K) | z), "log(q)")
  expect_identical(.response_name(.plain_response(log(q) ~ x | z)), ".sfa_y")
  expect_identical(.plain_response(y ~ x | z), y ~ x | z)
  d <- .add_response_column(log(q) ~ x, data.frame(q = c(1, exp(2)), x = 1:2))
  expect_equal(d[["log(q)"]], c(0, 2))
})

test_that("psfm() fits log(y) ~ x as the pre-transformed y ~ x", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 50)
  d$ey <- exp(d$y_ssfe)
  d$l_ey <- log(d$ey)
  d$eyt <- exp(d$y_tfe)
  d$l_eyt <- log(d$eyt)
  d$eyb <- exp(d$y_bc92)
  d$l_eyb <- log(d$eyb)
  expect_same_fit(
    psfm(log(ey) ~ x1 + x2, model_name = "PL80", data = d, individual = "name"),
    psfm(l_ey ~ x1 + x2, model_name = "PL80", data = d, individual = "name")
  )
  expect_same_fit(
    psfm(log(eyb) ~ x1 + x2, model_name = "BC92", data = d, individual = "name"),
    psfm(l_eyb ~ x1 + x2, model_name = "BC92", data = d, individual = "name")
  )
  ## The TFE family reads the response by name inside its own starting values.
  expect_same_fit(
    fit_tfe_quietly(psfm(log(eyt) ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d, individual = "name")),
    fit_tfe_quietly(psfm(l_eyt ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d, individual = "name"))
  )
})

test_that("psfm() on USUtilities: log() in the formula equals pre-logged columns", {
  skip_on_cran()
  data(USUtilities, package = "sfa", envir = environment())
  U <- as.data.frame(USUtilities)
  U <- transform(U, lq = log(q), lK = log(K), lL = log(L), lF = log(F))
  a <- psfm(log(q) ~ log(K) + log(L) + log(F),
    model_name = "PL80", data = U, individual = "firmID"
  )
  b <- psfm(lq ~ lK + lL + lF, model_name = "PL80", data = U, individual = "firmID")
  expect_same_fit(a, b)
  expect_lt(max(abs(a$coefficients)), 100)
})

test_that("sfm(), zsfm(), ttsfm() and lcsfm() fit a transformed response as transformed", {
  skip_on_cran()
  d <- cs_small()
  d$ey <- exp(d$y_pcs)
  d$l_ey <- log(d$ey)
  d$ez <- exp(d$y_pcs_z)
  d$l_ez <- log(d$ez)
  d$eyz <- exp(d$y_zisf)
  d$l_eyz <- log(d$eyz)
  d$eyt <- exp(d$y_ttne)
  d$l_eyt <- log(d$eyt)
  expect_same_fit(
    sfm(log(ey) ~ x1 + x2, model_name = "NHN", data = d),
    sfm(l_ey ~ x1 + x2, model_name = "NHN", data = d)
  )
  expect_same_fit(
    suppressWarnings(sfm(log(ez) ~ x1 + x2 | z, model_name = "NHN_Z", data = d)),
    suppressWarnings(sfm(l_ez ~ x1 + x2 | z, model_name = "NHN_Z", data = d))
  )
  expect_same_fit(
    suppressWarnings(zsfm(log(eyz) ~ x1 + x2, data = d)),
    suppressWarnings(zsfm(l_eyz ~ x1 + x2, data = d))
  )
  expect_same_fit(
    suppressWarnings(ttsfm(log(eyt) ~ x1 + x2, model_name = "TTNE", data = d)),
    suppressWarnings(ttsfm(l_eyt ~ x1 + x2, model_name = "TTNE", data = d))
  )
  expect_same_fit(
    suppressWarnings(lcsfm(log(ey) ~ x1 + x2, data = d)),
    suppressWarnings(lcsfm(l_ey ~ x1 + x2, data = d))
  )
})

test_that("residuals() of a log-response fit are on the log scale", {
  skip_on_cran()
  d <- cs_small()
  d$ey <- exp(d$y_pcs)
  d$l_ey <- log(d$ey)
  a <- sfm(log(ey) ~ x1 + x2, model_name = "NHN", data = d)
  b <- sfm(l_ey ~ x1 + x2, model_name = "NHN", data = d)
  expect_equal(residuals(a), residuals(b), tolerance = 1e-6)
})

test_that("gof_test() bootstraps a log-response fit on the fitted scale", {
  skip_on_cran()
  d <- cs_small()
  d$ey <- exp(d$y_pcs)
  d$l_ey <- log(d$ey)
  a <- sfm(log(ey) ~ x1 + x2, model_name = "NHN", data = d)
  b <- sfm(l_ey ~ x1 + x2, model_name = "NHN", data = d)
  expect_equal(gof_test(a, B = 5, seed = 1), gof_test(b, B = 5, seed = 1))
})
