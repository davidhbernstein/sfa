## The "sfareg" layout contract, asserted for EVERY model_name in the package.
##
## `out` is p x 3 -- one ROW per parameter, columns par/st_err/t-val -- and
## coefficients/std.errors/t.values are its columns. The source builds it
## 3 x p and stores t(out), so a branch that forgets the transpose returns an
## object on which the documented fit$out[, "par"] fails.
##
## That is not hypothetical: psfm("PL80_MVTN") shipped that way until 1.2.1.
## It survived because the contract was asserted per entry point, over a
## hand-written list of model_names, and PL80_MVTN was on none of them.
##
## So the lists below are DERIVED from each function's own formals, and the
## first test fails if a model_name exists with no spec here. Adding a model
## to the package therefore forces adding it to this file.

.contract <- function(fit, lbl) {
  expect_s3_class(fit, "sfareg")
  o <- fit$out
  expect_true(is.matrix(o), info = lbl)
  expect_identical(ncol(o), 3L, info = lbl)
  expect_identical(colnames(o), c("par", "st_err", "t-val"), info = lbl)
  expect_identical(rownames(o), names(fit$coefficients), info = lbl)
  expect_equal(unname(o[, "par"]), unname(fit$coefficients), info = lbl)
  expect_equal(unname(o[, "st_err"]), unname(fit$std.errors), info = lbl)
  expect_equal(unname(o[, "t-val"]), unname(fit$t.values), info = lbl)
  expect_identical(fit$model_name, sub("^[a-z]+ ", "", lbl))
}

## model_name -> the response column and formula that identifies it.
.sfm_specs <- list(
  NHN = "y_pcs ~ x1 + x2",        NHN_Z = "y_pcs_z ~ x1 + x2 | z",
  NE = "y_pcs_e ~ x1 + x2",       NE_Z  = "y_pcs_ez ~ x1 + x2 | z",
  NR = "y_pcs_r ~ x1 + x2",       THT   = "y_pcs_t ~ x1 + x2",
  NTN = "y_pcs_tn ~ x1 + x2",     NG    = "y_pcs_g ~ x1 + x2",
  NNAK = "y_pcs_nak ~ x1 + x2",   NU    = "y_pcs_u ~ x1 + x2",
  NGE = "y_pcs_ge ~ x1 + x2",     NLN   = "y_pcs_ln ~ x1 + x2",
  NW = "y_pcs_wb ~ x1 + x2",      tHN   = "y_pcs_thn ~ x1 + x2",
  TSL = "y_pcs_tsl ~ x1 + x2",    NGB2  = "y_pcs ~ x1 + x2",
  ## NB refuses to fit continuous data by design; its refusal is tested in
  ## test-nb.R. There is no fitted object here to check a layout on.
  NB = NA_character_
)

.psfm_specs <- list(
  TFE = "y_tfe ~ x1_w + x2_w",        TFE_WMLE = "y_tfe ~ x1_w + x2_w",
  SSFE = "y_ssfe ~ x1 + x2",          FD = "y_fd ~ x_fd | z_fd",
  TRE = "y_tre ~ x1 + x2",            GTRE = "y_gtre ~ x1 + x2",
  GTRE_FML = "y_gtre ~ x1 + x2",      PL80 = "y_ssfe ~ x1 + x2",
  PL80_MVTN = "y_pl_mvtn ~ x1 + x2",  BC92 = "y_bc92 ~ x1 + x2",
  K1990 = "y_bc92 ~ x1 + x2",         K1990modified = "y_bc92 ~ x1 + x2",
  TRE_Z = "y_tre_z ~ x1 + x2 | z_gtre", GTRE_Z = "y_gtre_z ~ x1 + x2 | z_gtre",
  GTRE_SEQ1 = "y_gtre ~ x1 + x2",     GTRE_SEQ2 = "y_gtre ~ x1 + x2",
  SSRE = "y_ssfe ~ x1 + x2",          SSCRE = "y_ssfe ~ x1 + x2",
  CSS = "y_ssfe ~ x1 + x2",           LS = "y_ssfe ~ x1 + x2",
  KSS = "y_ssfe ~ x1 + x2"
)

test_that("every model_name in the package has a layout spec in this file", {
  ## The guard that makes the rest of the file self-maintaining.
  expect_setequal(names(.sfm_specs), eval(formals(sfm)$model_name))
  expect_setequal(names(.psfm_specs), eval(formals(psfm)$model_name))
  expect_setequal(c("ZISF", "ZISF_Z"), eval(formals(zsfm)$model_name))
  expect_setequal(c("LCM", "LCM_Z", "LCM_CN"), eval(formals(lcsfm)$model_name))
  expect_setequal(c("TTNE", "TTHN", "TTNLS"), eval(formals(ttsfm)$model_name))
  expect_setequal(c("IVLIML", "IVCF", "C2SLS"), eval(formals(ivsfm)$model_name))
})

test_that("every sfm() model returns the documented p x 3 layout", {
  skip_on_cran()
  d <- as.data.frame(cs_small(N = 300))
  for (m in names(.sfm_specs)) {
    if (is.na(.sfm_specs[[m]])) next
    code <- sprintf("sfm(%s, model_name = '%s', data = d)", .sfm_specs[[m]], m)
    fit <- suppressWarnings(eval(parse(text = code)[[1]]))
    .contract(fit, paste("sfm", m))
  }
})

test_that("every psfm() model returns the documented p x 3 layout", {
  skip_on_cran()
  dp <- panel_small(t = 6, N = 50)
  for (m in names(.psfm_specs)) {
    extra <- if (identical(m, "GTRE")) ", estimator = 'sml'" else ""
    code <- sprintf("psfm(%s, model_name = '%s', data = dp, individual = 'name'%s)",
      .psfm_specs[[m]], m, extra)
    fit <- suppressWarnings(eval(parse(text = code)[[1]]))
    .contract(fit, paste("psfm", m))
  }
})

test_that("the remaining entry points return the documented p x 3 layout", {
  skip_on_cran()
  d <- as.data.frame(cs_small(N = 300))
  .contract(suppressWarnings(
    zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d)), "zsfm ZISF")
  .contract(suppressWarnings(
    zsfm(y_zisf_z ~ x1 + x2 | z, model_name = "ZISF_Z", data = d)), "zsfm ZISF_Z")
  .contract(suppressWarnings(
    lcsfm(y_pcs ~ x1 + x2, model_name = "LCM", data = d)), "lcsfm LCM")
  .contract(suppressWarnings(
    lcsfm(y_pcs ~ x1 + x2 | z, model_name = "LCM_Z", data = d)), "lcsfm LCM_Z")
  .contract(suppressWarnings(
    lcsfm(y_pcs ~ x1 + x2, model_name = "LCM_CN", data = d)), "lcsfm LCM_CN")
  for (m in c("TTNE", "TTHN", "TTNLS")) {
    .contract(suppressWarnings(
      ttsfm(y_ttne ~ x1 + x2, model_name = m, data = d)), paste("ttsfm", m))
  }
})
