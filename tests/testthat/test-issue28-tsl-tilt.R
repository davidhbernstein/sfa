## Issue #28: TSL's two exponential tilts, at a collapsing sigma_u.
##
## sfm()'s TSL likelihood built each tilt as two numbers,
##   A  = sig_v^2 / (2 sig_u^2) + eps / sig_u
##   l1 = log 2 + A + log Phi(-sig_v/sig_u - eps/sig_v),
## and with t = sig_v/sig_u, w = eps/sig_v that is
##   log 2 + t^2/2 + w t + log Phi(-(t + w)),
## whose t^2/2 + w t cancels EXACTLY against the same terms inside log Phi,
## because log phi(t + w) = -log sqrt(2 pi) - t^2/2 - t w - w^2/2. Added as two
## numbers they cancel catastrophically as sig_u -> 0. Measured on main before
## the fix: identical to the stable form through t = 1e4, off by 450 at t = 1e8,
## and returning +2.55e40 at sig_v = 1.126e20, sig_u = 1e-7 -- finite, hugely
## "better" than any real optimum, so the non-finite guard never saw it.

## The cancellation is an algebraic identity, checked here on its own so a
## failure separates "the algebra is wrong" from ".log_phi_tilt() is wrong".
test_that("the TSL tilts reduce to -eps^2 / (2 sig_v^2) exactly", {
  set.seed(11)
  eps <- rnorm(50, -0.5, 1)
  for (sig_v in c(0.2, 1, 3)) {
    for (sig_u in c(0.3, 1, 4)) {
      for (lam in c(0.2, 1, 5)) {
        A <- sig_v^2 / (2 * sig_u^2) + eps / sig_u
        B <- (1 + lam)^2 * sig_v^2 / (2 * sig_u^2) + eps * (1 + lam) / sig_u
        aa <- -sig_v / sig_u - eps / sig_v
        bb <- -sig_v * (1 + lam) / sig_u - eps / sig_v
        q <- eps^2 / (2 * sig_v^2)
        ## The residual is judged against the magnitude of the terms that
        ## cancelled, not against q: A - aa^2/2 is exactly -q in exact
        ## arithmetic, and what is left in floating point is the cancellation
        ## error this whole file is about. At lam = 5, sig_u = 0.3 the cancelled
        ## part of B is ~1.8e3, so ~1e-12 absolute IS machine precision here --
        ## which is the point: the error scales with the tilt, not with q.
        expect_lt(max(abs((A - aa^2 / 2) - (-q))),
          1e-13 * max(abs(A)) + 1e-14)
        expect_lt(max(abs((B - bb^2 / 2) - (-q))),
          1e-13 * max(abs(B)) + 1e-14)
      }
    }
  }
})

## Where nothing cancels, the shipped form must still agree with the naive one:
## this is what stops the rewrite from quietly changing the model.
test_that("the stable TSL tilts match the naive form at ordinary values", {
  set.seed(12)
  eps <- rnorm(200, -0.5, 1)
  naive <- function(eps, sig_v, sig_u, lam) {
    A <- sig_v^2 / (2 * sig_u^2) + eps / sig_u
    B <- (1 + lam)^2 * sig_v^2 / (2 * sig_u^2) + eps * (1 + lam) / sig_u
    l1 <- log(2) + A + stats::pnorm(-sig_v / sig_u - eps / sig_v, log.p = TRUE)
    l2 <- B + stats::pnorm(-sig_v * (1 + lam) / sig_u - eps / sig_v, log.p = TRUE)
    c(l1, l2)
  }
  stable <- function(eps, sig_v, sig_u, lam) {
    q <- eps^2 / (2 * sig_v^2)
    l1 <- log(2) - q + .log_phi_tilt(-sig_v / sig_u - eps / sig_v)
    l2 <- -q + .log_phi_tilt(-sig_v * (1 + lam) / sig_u - eps / sig_v)
    c(l1, l2)
  }
  for (sig_v in c(0.2, 0.5, 1, 2)) {
    for (sig_u in c(0.3, 1, 2)) {
      for (lam in c(0.2, 1, 5)) {
        expect_equal(stable(eps, sig_v, sig_u, lam),
          naive(eps, sig_v, sig_u, lam), tolerance = 1e-10,
          info = sprintf("sig_v=%g sig_u=%g lam=%g", sig_v, sig_u, lam))
      }
    }
  }
})

## The naive form's own failure, pinned so the test is about a real defect and
## not only about agreement. At t = sig_v/sig_u = 1e8 the two forms disagree by
## hundreds of log-likelihood units, and the naive one is the wrong one.
test_that("the naive TSL tilt really does lose the cancellation", {
  set.seed(13)
  eps <- rnorm(300, -0.5, 1)
  sig_v <- 1
  sig_u <- 1e-8
  lam <- 1
  A <- sig_v^2 / (2 * sig_u^2) + eps / sig_u
  aa <- -sig_v / sig_u - eps / sig_v
  naive_l1 <- log(2) + A + stats::pnorm(aa, log.p = TRUE)
  stable_l1 <- log(2) - eps^2 / (2 * sig_v^2) + .log_phi_tilt(aa)
  ## Both are finite -- which is exactly why the non-finite guard cannot help.
  expect_true(all(is.finite(naive_l1)))
  expect_true(all(is.finite(stable_l1)))
  expect_gt(max(abs(naive_l1 - stable_l1)), 1)
})

## End-to-end: the sample from issue #28. TSL on y_pcs_thn is a misspecified
## fit with no inefficiency signal, so sigma_u collapses and the correct answer
## is the OLS log-likelihood. Before the fix sfm() reported +1.70141183e+40 at
## sig_v = 1.1259e20, sig_u = 1e-7, lambda = 6.528e6, with efficiency
## identically 1.
test_that("TSL does not return a non-physical fit on the issue #28 sample", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 1002715, cons = .5, beta1 = .5, beta2 = .5,
    sig_u = 1, sig_v = 1, mu = .5, a = 5)
  d$yy <- d$y_pcs_thn
  ll_ols <- as.numeric(stats::logLik(stats::lm(yy ~ x1 + x2, data = d)))
  f <- sfm(yy ~ x1 + x2, data = d, model_name = "TSL")
  ll <- as.numeric(logLik(f))

  expect_true(is.finite(ll))
  ## The defect was a spurious optimum 1e40 units "better" than any real one.
  expect_lt(abs(ll), 1e4)
  ## On this sample TSL's profile likelihood rises monotonically to OLS as
  ## sigma_u -> 0 (profile - OLS: -4.18 at sigma_u = 1, -3.6e-3 at 0.1,
  ## -5.4e-6 at 0.01, -5.6e-9 at 1e-3), so OLS is the supremum, reached only
  ## on the boundary, and the fit cannot exceed it. How close the optimizer
  ## gets along that flat approach is platform-dependent -- 9e-4 short on
  ## CI's Linux and Windows runners -- so only the bound is asserted.
  ##
  ## THE PRECONDITION IS WRONG SKEW, AND DROPPING IT WOULD BE WRONG. OLS is the
  ## supremum only because these residuals are skewed the WRONG way (third
  ## central moment +0.541 here), so there is no interior maximum. On a
  ## correctly skewed sample TSL's interior MLE legitimately BEATS OLS, and by
  ## a lot: on y_pcs_tsl, TSL's own DGP, all 24 fits measured exceed OLS, by up
  ## to 21.3 units. Across 432 fits spanning 18 response columns, 262 exceed
  ## OLS by more than 1 unit. So `ll <= ll_ols` is NOT a general property of
  ## the model and must never be asserted unconditionally -- it holds here
  ## because of the skew check below, which is what makes this sample a
  ## boundary case.
  ## Wrong skew is a property of the OLS RESIDUALS, not of the raw response --
  ## the package's own .warn_wrong_skew_boundary() reports the residual moment,
  ## and on this sample that is +0.541, the figure it prints. The two are not
  ## interchangeable: over the 432 (seed, column) pairs used above, the sign of
  ## m3(y) and the sign of m3(resid) disagree on 28, i.e. 6.5%. Both happen to
  ## agree on this seed, so testing the raw response would pass here and be
  ## wrong for anyone who adds a sample.
  expect_gt(mean(resid(stats::lm(yy ~ x1 + x2, data = d))^3), 0)
  expect_lte(ll, ll_ols + 1e-6)
  ## What generalises instead is the CONJUNCTION: sigma_u pinned at its floor
  ## AND a likelihood above OLS. At the floor the correct likelihood is within
  ## ~1e-11 of OLS whatever the skew, so that pair is unreachable by any
  ## correct fit. Of the 262 fits above that legitimately exceed OLS, ZERO sit
  ## on the floor; of the 432, the conjunction flags 17 before this fix and 0
  ## after, catching every non-physical fit and every quiet 1-11 unit error.
  expect_false(f$out["sigu", "par"] <= 1.0000001e-07 && ll > ll_ols + 1e-6)
  ## sig_v ran to 1.1259e20 before the fix.
  expect_lt(f$out["sigv", "par"], 1e3)
  expect_lt(f$out["lambda", "par"], 1e4)
})
