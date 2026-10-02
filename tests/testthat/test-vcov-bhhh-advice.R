## vcov(type = "sandwich"/"clustered") advises type = "bhhh" when its own
## bread is undefined. That advice has to be true.
##
## Both failure paths used to recommend BHHH unconditionally, reasoning that it
## needs no Hessian. True, and not sufficient: BHHH still inverts
## crossprod(G), and the fits that reach these errors are often flat in a
## parameter, which kills the outer product as well. On the fit below the old
## message said `type = "bhhh"` "is defined here" while `vcov(type = "bhhh")`
## raised "the outer product of gradients is singular" -- so a user following
## the advice hit a second error.
##
## The trigger is a real degenerate fit rather than a contrived matrix: NTN on
## these samples has no interior maximum in lambda. The profile log-likelihood
## rises monotonically to a finite limit as lambda -> Inf, so the likelihood
## really is flat in lambda at the reported fit and the Hessian dies with it.
##
## WHICH of them also makes the SANDWICH bread singular is knife-edge, because
## the fit sits on that flat ridge. Issue #55 demonstrated this: a correct
## change to the NTN likelihood slid the rand = 2 fit from lambda = 9.7e5 to
## 4.4e6 -- a log-likelihood move of 9e-5 along the ridge -- and the bread
## became invertible again while solve(hessian) still failed and BHHH still
## failed. The invariant was never in danger; the trigger simply stopped
## triggering. So the trigger is SEARCHED FOR over several samples rather than
## hard-coded, and the test fails if none of them exercises it.

test_that("the bhhh recommendation agrees with whether bhhh works", {
  skip_on_cran()
  cands <- lapply(c(2, 4, 5, 12), function(sd) {
    d <- data_gen_cs(N = 150, rand = sd, cons = .5, beta1 = .5, beta2 = .5,
      sig_u = 1, sig_v = .5, mu = .5, a = 5)
    list(d = d, f = suppressWarnings(sfm(y_pcs_tn ~ x1 + x2, data = d,
      model_name = "NTN", keep_objective = TRUE)))
  })

  ## Every candidate must be the kind of fit this test is about: lambda has run
  ## away and no standard errors exist.
  for (cd in cands) {
    expect_gt(cd$f$out["lambda", "par"], 1e4)
    expect_true(all(is.na(cd$f$out[, "st_err"])))
  }

  exercised <- 0L
  for (cd in cands) {
    f <- cd$f
    bhhh <- tryCatch(vcov(f, type = "bhhh"), error = function(e) e)
    bhhh_works <- !inherits(bhhh, "error")

    for (ty in c("sandwich", "clustered")) {
      e <- tryCatch(vcov(f, type = ty, cluster = seq_len(nrow(cd$d))),
        error = function(e) e)
      ## The bread is defined on this fit, so there is no advice to check.
      if (!inherits(e, "error")) next
      exercised <- exercised + 1L
      msg <- conditionMessage(e)

      ## THE INVARIANT: the message may claim BHHH is available only when it is.
      ## Asserted as an agreement between two observed behaviours, so it cannot
      ## pass by the wording happening to change.
      claims_available <- grepl("is defined here", msg, fixed = TRUE)
      expect_equal(claims_available, bhhh_works,
        info = sprintf("type = %s: message claims bhhh available = %s, actual = %s",
          ty, claims_available, bhhh_works))

      ## When BHHH is genuinely unavailable the message must say so rather than
      ## stay silent about it.
      if (!bhhh_works) {
        expect_match(msg, "ALSO undefined", fixed = TRUE)
      }
      ## And either way it still has to name the parameter that died (gap A55).
      expect_match(msg, "lambda")
    }
  }

  ## A search that found nothing would make every assertion above vacuous.
  expect_gt(exercised, 0L)
})

test_that("a healthy fit still returns all three covariances", {
  skip_on_cran()
  ## Guards the other direction: the availability probe added for the above
  ## must not disturb the normal path.
  d <- data_gen_cs(N = 200, rand = 11, cons = .5, beta1 = .5, beta2 = .5,
    sig_u = 1, sig_v = .5, mu = .5, a = 5)
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = d,
    model_name = "NHN", keep_objective = TRUE))
  p <- length(coef(f))
  for (ty in c("hessian", "bhhh", "sandwich")) {
    V <- vcov(f, type = ty)
    expect_equal(dim(V), c(p, p), info = ty)
    expect_true(all(is.finite(diag(V))), info = ty)
  }
  V <- vcov(f, type = "clustered", cluster = seq_len(nrow(d)))
  expect_equal(dim(V), c(p, p))
})
