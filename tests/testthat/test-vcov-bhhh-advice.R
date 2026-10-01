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
## this sample has no interior maximum in lambda. Its profile log-likelihood
## rises monotonically to a finite limit as lambda -> Inf (-188.1003 at
## lambda = 9.7e5, flat at -188.1001876 from 1e9 to 1e15), so the likelihood
## really is flat in lambda at the reported fit, and the score column for
## lambda is 1.24e-10 against 1.15e+07 for the largest. crossprod(G) has
## condition number 6.8e32.

test_that("the bhhh recommendation agrees with whether bhhh works", {
  skip_on_cran()
  d <- data_gen_cs(N = 150, rand = 2, cons = .5, beta1 = .5, beta2 = .5,
    sig_u = 1, sig_v = .5, mu = .5, a = 5)
  f <- suppressWarnings(sfm(y_pcs_tn ~ x1 + x2, data = d,
    model_name = "NTN", keep_objective = TRUE))

  ## The fit this test needs: lambda has run away and no standard errors exist.
  expect_gt(f$out["lambda", "par"], 1e4)
  expect_true(all(is.na(f$out[, "st_err"])))

  bhhh <- tryCatch(vcov(f, type = "bhhh"), error = function(e) e)
  bhhh_works <- !inherits(bhhh, "error")

  for (ty in c("sandwich", "clustered")) {
    e <- tryCatch(vcov(f, type = ty, cluster = seq_len(nrow(d))),
      error = function(e) e)
    expect_s3_class(e, "error")
    msg <- conditionMessage(e)

    ## THE INVARIANT: the message may claim BHHH is available only when it is.
    ## Asserted as an agreement between two observed behaviours, so it cannot
    ## pass by the wording happening to change.
    claims_available <- grepl("bhhh\\\\\"? needs no Hessian and is defined here",
      msg) || grepl("is defined here", msg)
    expect_equal(claims_available, bhhh_works,
      info = sprintf("type = %s: message claims bhhh available = %s, actual = %s",
        ty, claims_available, bhhh_works))

    ## On this fit BHHH is genuinely unavailable, so the message must say so
    ## rather than stay silent about it.
    if (!bhhh_works) {
      expect_match(msg, "ALSO undefined", fixed = TRUE)
    }
    ## And either way it still has to name the parameter that died (gap A55).
    expect_match(msg, "lambda")
  }
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
