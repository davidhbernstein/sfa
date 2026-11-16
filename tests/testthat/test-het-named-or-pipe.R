## A determinant block may be written as a pipe segment or as a named
## argument (1.3.0). The two are the same thing: the named form is spliced
## into its pipe position before anything else runs.
##
## THE FIRST THING ASSERTED IS THAT PIPES DID NOT CHANGE. The alias layer is
## only reachable when a named argument is actually supplied, and the tests
## below pin that both as an identity on the formula and as bit-identical
## fits, because "nothing breaks" is the property that matters most here.

cs <- function(n = 120, seed = 3) {
  data_gen_cs(N = n, rand = seed, sig_u = 1, sig_v = .3, cons = .5,
              beta1 = .5, beta2 = .5, a = 5, mu = .5)
}

## ---- pipes are untouched ---------------------------------------------------

test_that("with no named argument the formula is returned UNCHANGED", {
  ## Not merely equal -- the identical object, so a pipe call provably does not
  ## enter the new code path.
  f <- y ~ x1 + x2 | z
  for (m in c("NHN_Z", "TRE_Z", "GTRE_Z", "ZISF_Z", "LCM_Z", "TTNE", "FD")) {
    expect_identical(
      .sfa_splice_het_formula(f, m, list(uhet = NULL, whet = NULL,
                                         hhet = NULL, ghet = NULL,
                                         classhet = NULL, scaling = NULL),
                              "test"),
      f, info = m)
  }
})

test_that("every pipe form still fits", {
  skip_on_cran()
  d <- cs()
  expect_s3_class(ttsfm(y_ttne ~ x1 + x2 | z | zp, data = d, model_name = "TTNE"), "sfareg")
  expect_s3_class(ttsfm(y_ttne ~ x1 + x2 | z, data = d, model_name = "TTNE"), "sfareg")
  expect_s3_class(suppressWarnings(
    zsfm(y_zisf_z ~ x1 + x2 | z, data = d, model_name = "ZISF_Z")), "sfareg")
  ## and a model with NO determinant block is unaffected either way
  expect_s3_class(ttsfm(y_ttne ~ x1 + x2, data = d, model_name = "TTNE"), "sfareg")
})

## ---- the named form is the same model --------------------------------------

test_that("named and pipe give bit-identical fits", {
  skip_on_cran()
  d <- cs()
  pairs <- list(
    list("TTNE uhet+whet",
         function() ttsfm(y_ttne ~ x1 + x2 | z | zp, data = d, model_name = "TTNE"),
         function() ttsfm(y_ttne ~ x1 + x2, data = d, model_name = "TTNE",
                          uhet = ~z, whet = ~zp)),
    list("TTNE uhet only",
         function() ttsfm(y_ttne ~ x1 + x2 | z, data = d, model_name = "TTNE"),
         function() ttsfm(y_ttne ~ x1 + x2, data = d, model_name = "TTNE", uhet = ~z)),
    list("TTHN uhet+whet",
         function() ttsfm(y_tthn ~ x1 + x2 | z | zp, data = d, model_name = "TTHN"),
         function() ttsfm(y_tthn ~ x1 + x2, data = d, model_name = "TTHN",
                          uhet = ~z, whet = ~zp)),
    list("ZISF_Z ghet",
         function() zsfm(y_zisf_z ~ x1 + x2 | z, data = d, model_name = "ZISF_Z"),
         function() zsfm(y_zisf_z ~ x1 + x2, data = d, model_name = "ZISF_Z", ghet = ~z))
  )
  for (p in pairs) {
    a <- suppressWarnings(p[[2]]()); b <- suppressWarnings(p[[3]]())
    ## tolerance = 0: the two routes are one code path, so anything other than
    ## exact equality means the splice changed the model.
    expect_equal(a$out[, "par"], b$out[, "par"], tolerance = 0, info = p[[1]])
    expect_identical(deparse(a$formula), deparse(b$formula), info = p[[1]])
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 0,
                 info = p[[1]])
  }
})

test_that("a second-slot argument on its own means `| 1 | zp`", {
  ## An omitted segment is `1`, which is what the pipe syntax already means,
  ## so whet without uhet must not be an error.
  f <- .sfa_splice_het_formula(y ~ x1 + x2, "TTNE",
                               list(uhet = NULL, whet = ~zp), "ttsfm")
  expect_identical(deparse(f), deparse(y ~ x1 + x2 | 1 | zp))
})

test_that("a transformed response survives the splice", {
  f <- .sfa_splice_het_formula(log(q) ~ x1, "NHN_Z", list(uhet = ~z), "sfm")
  expect_identical(deparse(f), deparse(log(q) ~ x1 | z))
})

## ---- refusals --------------------------------------------------------------

test_that("giving a block both ways is refused, not silently resolved", {
  expect_error(
    .sfa_splice_het_formula(y ~ x | z, "NHN_Z", list(uhet = ~z2), "sfm"),
    "given twice")
  expect_error(
    .sfa_splice_het_formula(y ~ x | z | zp, "GTRE_Z", list(hhet = ~z2), "psfm"),
    "given twice")
})

test_that("a named argument is refused where its pipe segment is", {
  ## The alias layer must not add a determinant block to a model that has
  ## none -- that would be a new model, not a new spelling.
  e <- tryCatch(.sfa_splice_het_formula(y ~ x, "GTRE", list(uhet = ~z), "psfm"),
                error = function(z) conditionMessage(z))
  expect_match(e, "not .GTRE.")
  expect_match(e, "NHN_Z")          # names the models that do have it
  expect_error(.sfa_splice_het_formula(y ~ x, "TTNE", list(hhet = ~z), "ttsfm"),
               "GTRE_Z")
})

test_that("a non-formula determinant is refused with the right shape shown", {
  expect_error(.sfa_splice_het_formula(y ~ x, "NHN_Z", list(uhet = "z"), "sfm"),
               "one-sided formula")
  expect_error(.sfa_splice_het_formula(y ~ x, "NHN_Z", list(uhet = y ~ z), "sfm"),
               "one-sided formula")
})

test_that("the entry points expose the arguments the slot map says they do", {
  expect_true(all(c("uhet", "hhet", "scaling") %in% names(formals(psfm))))
  expect_true(all(c("uhet", "whet") %in% names(formals(ttsfm))))
  expect_true("ghet" %in% names(formals(zsfm)))
  expect_true("classhet" %in% names(formals(lcsfm)))
  ## and they are optional, so every existing call is unaffected
  for (a in c("uhet", "hhet", "scaling")) expect_null(eval(formals(psfm)[[a]]))
  for (a in c("uhet", "whet")) expect_null(eval(formals(ttsfm)[[a]]))
  expect_null(eval(formals(zsfm)$ghet))
  expect_null(eval(formals(lcsfm)$classhet))
})

## ---- the stored formula, and the stored call --------------------------------

test_that("$formula carries the determinants and $call records what was typed", {
  skip_on_cran()
  d <- cs()
  f <- ttsfm(y_ttne ~ x1 + x2, data = d, model_name = "TTNE", uhet = ~z)
  g <- ttsfm(y_ttne ~ x1 + x2 | z, data = d, model_name = "TTNE")
  ## $formula is the SPLICED one, because that is what every consumer reads --
  ## psfm_bootstrap(), .sfa_xb(), residuals(), lcsfm_homogeneity(). If it were
  ## the unspliced one those would silently lose the determinant block.
  ##
  ## Asserted against the PIPE ROUTE rather than against a literal, because
  ## ttsfm() normalises a two-part formula to three parts itself -- both routes
  ## store `y_ttne ~ x1 + x2 | z | 1`. Hard-coding the caller's own two-part
  ## formula here failed for that reason, and comparing the two routes is the
  ## property that actually matters anyway: it cannot go stale if the
  ## normalisation changes.
  expect_identical(deparse(f$formula), deparse(g$formula))
  expect_true(grepl("\\| *z", paste(deparse(f$formula), collapse = "")))
  ## $call keeps the caller's own syntax, which is the honest record of what
  ## was asked for and is what print()/spec_test() label the fit with.
  expect_true("uhet" %in% names(as.list(f$call)))
  expect_identical(deparse(f$call$formula), deparse(y_ttne ~ x1 + x2))
})

## ---- psfm and lcsfm end to end ---------------------------------------------

test_that("psfm and lcsfm named routes also give bit-identical fits", {
  skip_on_cran()
  dp <- data_gen_p(N = 30, t = 5, rand = 2, sig_u = 1, sig_v = .5, sig_r = .5,
                   sig_h = .5, cons = .5, beta1 = .5, beta2 = .5)
  a <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, data = dp,
                             model_name = "TRE_Z", individual = "name",
                             halton_num = 25))
  b <- suppressWarnings(psfm(y_tre_z ~ x1 + x2, data = dp, model_name = "TRE_Z",
                             individual = "name", halton_num = 25,
                             uhet = ~z_gtre))
  expect_equal(a$out[, "par"], b$out[, "par"], tolerance = 0, info = "psfm TRE_Z uhet")
  expect_identical(deparse(a$formula), deparse(b$formula))

  dc <- cs()
  p <- suppressWarnings(lcsfm(y_pcs ~ x1 + x2 | z, data = dc, model_name = "LCM_Z"))
  q <- suppressWarnings(lcsfm(y_pcs ~ x1 + x2, data = dc, model_name = "LCM_Z",
                              classhet = ~z))
  expect_equal(p$out[, "par"], q$out[, "par"], tolerance = 0, info = "lcsfm LCM_Z classhet")
  expect_identical(deparse(p$formula), deparse(q$formula))
})
