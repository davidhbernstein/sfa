## Reported parameter names have to identify which equation a parameter
## belongs to (gaps A80 and A84).
##
## A80: the variance-determinant design matrix carries its own intercept, so a
## homoskedastic `| z` segment gave two rows both called "(Intercept)" -- one
## from the frontier, one from the sigma_u equation. Nothing in print(),
## summary() or a texreg table distinguished them, and any name-based lookup on
## $out self-matched silently.
##
## A84: the prefixed names 1.3.0 introduced were not recognised as variance
## parameters at all, so a texreg table of an `sfm(vhet = ~z)` fit put them in
## the frontier body and include.scales = FALSE did not drop them.
##
## The assertions are the INVARIANTS -- no duplicate names, and every
## determinant row classified as a scale -- not the specific strings, so a
## later change of convention does not have to come here to be re-blessed.

zdata <- function(seed = 6, n = 300) {
  data_gen_cs(N = n, rand = seed, sig_u = 1, sig_v = .3, cons = .5,
              beta1 = .5, beta2 = .5, a = 5, mu = .5)
}

## Every entry point that takes a variance-determinant pipe segment and is
## cheap enough to fit here. The panel _Z models are covered by psfm.R's own
## de-duplication (R/psfm.R:1321, :1775) and cost minutes to fit, so they are
## checked by that code's own tests rather than refitted here.
pipe_fits <- list(
  NHN_Z  = function(d) sfm(y_pcs_z ~ x1 + x2 | z, data = d, model_name = "NHN_Z"),
  NE_Z   = function(d) sfm(y_pcs_ez ~ x1 + x2 | z, data = d, model_name = "NE_Z"),
  ZISF_Z = function(d) zsfm(y_zisf_z ~ x1 + x2 | z, data = d, model_name = "ZISF_Z")
)

test_that("no reported parameter name is duplicated on a pipe model", {
  skip_on_cran()
  d <- zdata()
  for (nm in names(pipe_fits)) {
    f <- suppressWarnings(pipe_fits[[nm]](d))
    rn <- rownames(f$out)
    expect_equal(anyDuplicated(rn), 0L,
                 info = paste(nm, "rows:", paste(rn, collapse = " ")))
    ## A duplicated name also makes a name-based lookup unsafe, so assert the
    ## thing that actually depends on it: vcov() must line up with $out.
    expect_identical(rownames(vcov(f)), rn, info = nm)
    expect_identical(names(coef(f)), rn, info = nm)
  }
})

test_that("the determinant rows are recognised as variance parameters", {
  skip_on_cran()
  d <- zdata()
  for (nm in names(pipe_fits)) {
    f <- suppressWarnings(pipe_fits[[nm]](d))
    rn <- rownames(f$out)
    det_rows <- grep("^(Zv|Zu|Zw|Zmu|Zgamma|scale|shape)\\.", rn, value = TRUE)
    ## There IS a determinant block on each of these: an empty set here would
    ## make the expectation below vacuously true.
    expect_gt(length(det_rows), 0)
    expect_true(all(.sfa_is_scale(det_rows)), info = nm)
    ## and the frontier coefficients must NOT be swept up with them
    expect_false(.sfa_is_scale("(Intercept)"))
    expect_false(.sfa_is_scale("x1"))
  }
})

test_that("include.scales = FALSE drops the heteroskedasticity coefficients", {
  skip_on_cran()
  skip_if_no_texreg <- function() testthat::skip_if_not_installed("texreg")
  skip_if_no_texreg()
  set.seed(9)
  n <- 300
  zv <- runif(n, -1, 1); x1 <- rnorm(n); x2 <- rnorm(n)
  y <- .5 + .5 * x1 + .5 * x2 + rnorm(n, 0, exp(-1 + .5 * zv)) - abs(rnorm(n, 0, 1))
  h <- suppressWarnings(sfm(y ~ x1 + x2, vhet = ~zv,
                            data = data.frame(y, x1, x2, zv), model_name = "NHN"))
  ## The fit really does carry a Zv block, or the test proves nothing.
  expect_true(any(grepl("^Zv\\.", rownames(h$out))))

  body_with    <- extract.sfareg(h)@coef.names
  body_without <- extract.sfareg(h, include.scales = FALSE)@coef.names
  expect_true(any(grepl("^Z", body_with)))
  expect_false(any(grepl("^Z", body_without)))
  ## and the frontier survives either way
  expect_true(all(c("x1", "x2") %in% body_without))
})

test_that(".sfa_is_scale knows every prefix the package emits", {
  ## Derived from the sites that build these names: het.R (Zv./Zu./Zmu./scale.),
  ## sfm.R:345 (shape.), ttsfm()'s .tt_report (Zu./Zw.) and
  ## starting.values.R (Zu./Zgamma.).
  for (p in c("Zv", "Zu", "Zw", "Zmu", "Zgamma", "scale", "shape")) {
    expect_true(.sfa_is_scale(paste0(p, ".(Intercept)")), info = p)
    expect_true(.sfa_is_scale(paste0(p, ".some_covariate")), info = p)
  }
  ## The prefix test requires the dot, so an ordinary covariate whose name
  ## merely starts with one of these words is not swept up by it.
  expect_false(.sfa_is_scale("Zugzwang"))
  expect_false(.sfa_is_scale("scaled_input"))

  ## KNOWN, PRE-EXISTING false positive, asserted so it is visible rather than
  ## discovered again: the ORIGINAL unanchored pattern
  ## grepl("^(sigma|sig|...|shape)", nm) matches any name merely STARTING with
  ## one of those words, so a covariate called "shapely" or "sigmoid" is still
  ## classified as a scale parameter. Not introduced by the A84 fix and not
  ## widened by it; left alone here because anchoring those patterns risks
  ## reclassifying real parameters such as "sigmaSq" and "logit0", and the
  ## consequence is cosmetic (a row moves out of the frontier body).
  expect_true(.sfa_is_scale("shapely"))
  expect_true(.sfa_is_scale("sigmoid"))
})
