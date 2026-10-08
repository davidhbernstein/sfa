## A60. sigma_u on the zero boundary reads in OPPOSITE directions depending on
## the residual skew, and both the fit-time warning and the printed NOTE used to
## treat only one of them. Under WRONG skew the boundary is the correct MLE
## (Olson, Schmidt and Waldman 1980 Type I; Waldman 1982). Under CORRECT skew
## the data say inefficiency is present, so the boundary should not be the
## answer and the fit is suspect. The old code warned only in the first case and
## printed the Waldman exculpation in BOTH.
##
## These assertions drive the two reporting paths directly rather than through a
## fit, because the models that reach this state are on a flat ridge where the
## optimizer's stopping point is not a stable property across platforms -- the
## lesson from the vcov fixture and from issue #55's CI failures.

test_that("the printed NOTE does not cite Waldman under correct skew", {
  mk <- function(wrong) {
    structure(list(sigma_u_at_bound = TRUE, wrong_skew = wrong),
      class = "sfareg")
  }

  wrong <- paste(capture.output(sfa:::.sfa_report_boundary(mk(TRUE))),
    collapse = " ")
  right <- paste(capture.output(sfa:::.sfa_report_boundary(mk(FALSE))),
    collapse = " ")

  ## Both still say what happened.
  expect_match(wrong, "sigma_u is on the zero boundary")
  expect_match(right, "sigma_u is on the zero boundary")

  ## Wrong skew keeps the Waldman reading, which is correct there.
  expect_match(wrong, "wrong\\s+skew")
  expect_match(wrong, "correct MLE")
  expect_match(wrong, "Waldman 1982")

  ## Correct skew must NOT be told the boundary is the correct MLE. This is the
  ## assertion that fails on main: the old text printed the "Under wrong
  ## skewness this is the correct MLE, not a failure" sentence unconditionally,
  ## having only dropped the clause naming the skew.
  expect_false(grepl("correct MLE", right, fixed = TRUE))
  expect_false(grepl("not a failure", right, fixed = TRUE))

  ## And it must say the opposite thing: the fit is suspect.
  expect_match(right, "evidence that inefficiency")
  expect_match(right, "may not be a maximum")
  expect_match(right, "does NOT apply")
})

test_that("the boundary NOTE is silent when sigma_u is interior", {
  for (wrong in list(TRUE, FALSE, NA)) {
    o <- structure(list(sigma_u_at_bound = FALSE, wrong_skew = wrong),
      class = "sfareg")
    expect_identical(capture.output(sfa:::.sfa_report_boundary(o)), character(0))
  }
  ## NA is the "could not be computed" case from .wrong_skew_boundary(); it must
  ## not take either branch.
  o <- structure(list(sigma_u_at_bound = NA, wrong_skew = FALSE),
    class = "sfareg")
  expect_identical(capture.output(sfa:::.sfa_report_boundary(o)), character(0))
})

test_that(".warn_boundary_correct_skew() reports the moment and the code", {
  ws <- list(wrong_skew = FALSE, at_bound = TRUE, m3 = -1585.5)

  w <- tryCatch(sfa:::.warn_boundary_correct_skew(ws, "NW", "sigu", conv = 52L),
    warning = function(e) conditionMessage(e))
  expect_type(w, "character")
  expect_match(w, "sfm(model_name = \"NW\")", fixed = TRUE)
  expect_match(w, "sigu has collapsed to the boundary", fixed = TRUE)
  ## The residual third moment is the evidence, so it has to be in the message.
  expect_match(w, as.character(signif(ws$m3, 3)), fixed = TRUE)
  ## A non-zero convergence code is the actionable part of A60's three NW fits.
  expect_match(w, "convergence code 52", fixed = TRUE)
  expect_match(w, "did not converge", fixed = TRUE)
  ## It must not repeat the wrong-skew exculpation as if it applied, and must
  ## say so explicitly. "does NOT apply" is the PRINT path's wording; the
  ## warning denies the same reading in its own words, so assert the warning's.
  expect_false(grepl("correct maximum likelihood estimate", w, fixed = TRUE))
  expect_match(w, "Do NOT read it as no evidence of inefficiency", fixed = TRUE)
  expect_match(w, "applies only under wrong skew", fixed = TRUE)

  ## conv 0 and conv NA both drop the convergence clause rather than printing
  ## "code 0" or "code NA".
  for (cv in list(0L, NA_integer_)) {
    w0 <- tryCatch(sfa:::.warn_boundary_correct_skew(ws, "NW", "sigu", conv = cv),
      warning = function(e) conditionMessage(e))
    expect_false(grepl("convergence code", w0, fixed = TRUE))
    expect_match(w0, "may not be a maximum", fixed = TRUE)
  }
  ## The default is NA, i.e. no claim about convergence.
  wd <- tryCatch(sfa:::.warn_boundary_correct_skew(ws, "NW", "sigu"),
    warning = function(e) conditionMessage(e))
  expect_false(grepl("convergence code", wd, fixed = TRUE))
})

test_that("the two boundary warnings stay distinguishable", {
  ## A caller filtering on one must not catch the other.
  wsw <- list(wrong_skew = TRUE, at_bound = TRUE, m3 = 12.3)
  wsr <- list(wrong_skew = FALSE, at_bound = TRUE, m3 = -12.3)

  a <- tryCatch(sfa:::.warn_wrong_skew_boundary(wsw, "NE", "sigu"),
    warning = function(e) conditionMessage(e))
  b <- tryCatch(sfa:::.warn_boundary_correct_skew(wsr, "NE", "sigu"),
    warning = function(e) conditionMessage(e))

  expect_false(identical(a, b))
  ## The old message's defining claim is that the boundary IS the estimate.
  expect_match(a, "boundary IS the maximum likelihood estimate", fixed = TRUE)
  expect_false(grepl("boundary IS the maximum likelihood estimate", b,
    fixed = TRUE))
  ## The new one's defining claim is the opposite.
  expect_match(b, "unlikely to be the maximum likelihood estimate", fixed = TRUE)
  expect_false(grepl("unlikely to be the maximum likelihood estimate", a,
    fixed = TRUE))
})

test_that("tHN keeps its own reading of a collapsed sigma_u", {
  ## tHN warns separately that heavy-tailed noise absorbing the one-sided
  ## component is a known property of the model, so calling the same fit
  ## "suspect" would contradict it on the same console. tHN's row names put it
  ## through the generic block, so it has to be excluded explicitly.
  thn_both <- structure(list(sigma_u_at_bound = TRUE, wrong_skew = FALSE,
    model_name = "tHN", thn_sigma_u_at_bound = TRUE), class = "sfareg")
  o <- paste(capture.output(sfa:::.sfa_report_boundary(thn_both)),
    collapse = " ")
  expect_match(o, "known property of tHN")
  expect_false(grepl("may not be a maximum", o, fixed = TRUE))
  expect_false(grepl("does NOT apply", o, fixed = TRUE))
  ## And exactly one NOTE, not the tHN one plus a contradicting one.
  expect_equal(lengths(regmatches(o, gregexpr("NOTE:", o, fixed = TRUE)))[[1]], 1L)

  ## The two thresholds differ (tHN uses sig_u/sqrt(sig_u^2+sig_v^2) < 1e-3,
  ## the generic one sigu/sd(resid) < 1e-2), so the generic flag can fire alone.
  ## That must still report the boundary, with tHN's reading and not the
  ## suspect one.
  thn_gen <- structure(list(sigma_u_at_bound = TRUE, wrong_skew = FALSE,
    model_name = "tHN", thn_sigma_u_at_bound = FALSE), class = "sfareg")
  g <- paste(capture.output(sfa:::.sfa_report_boundary(thn_gen)), collapse = " ")
  expect_match(g, "sigma_u is on the zero boundary")
  expect_match(g, "known property")
  expect_false(grepl("may not be a maximum", g, fixed = TRUE))
  expect_false(grepl("correct MLE", g, fixed = TRUE))

  ## Under WRONG skew tHN still gets the Waldman reading, which is unchanged.
  thn_wrong <- structure(list(sigma_u_at_bound = TRUE, wrong_skew = TRUE,
    model_name = "tHN", thn_sigma_u_at_bound = FALSE), class = "sfareg")
  w <- paste(capture.output(sfa:::.sfa_report_boundary(thn_wrong)),
    collapse = " ")
  expect_match(w, "correct MLE")
})

test_that("a model with no model_name still gets the generic reading", {
  ## Defensive: identical(NULL, "tHN") is FALSE, so an object missing
  ## model_name must not silently fall into the tHN branch.
  o <- structure(list(sigma_u_at_bound = TRUE, wrong_skew = FALSE),
    class = "sfareg")
  out <- paste(capture.output(sfa:::.sfa_report_boundary(o)), collapse = " ")
  expect_match(out, "may not be a maximum")
})
