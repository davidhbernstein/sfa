## Both routes to the four-component model must report the SAME estimand.
##
## Until 1.3.0, estimator = "fiml" predicted transient efficiency with
## .te_battese_coelli() -- the two-component normal/half-normal posterior --
## applied to the raw residual, which still carries the individual effect and
## the persistent inefficiency. Neither sigma_r nor sigma_h entered it, so
## persistent inefficiency was charged to the transient term. Both routes now
## use .gtre_csn_te(), the closed-skew-normal posterior of Colombi, Kumbhakar,
## Martini and Vittadini (2014).
##
## These tests assert the invariants the fix guarantees, not where an
## optimizer lands.

test_that(".gtre_csn_te() reduces to the two-component posterior when the persistent components vanish", {
  ## The arithmetic is pinned by a LIMIT, not by a fit. With sigma_h and
  ## sigma_r driven to zero the four-component model collapses to
  ## noise-plus-transient-inefficiency, where the closed form is exactly
  ## .te_battese_coelli() on the residual. Nothing here is simulated.
  t_i <- 4L; N <- 3L
  e <- list(c(-0.9, -0.2, -1.6, -0.05),
            c(-0.3, -1.1, -0.4, -2.0),
            c(-0.01, -3.0, -0.7, -0.25))
  sig_u <- 0.9; sig_v <- 0.35
  tiny  <- 1e-6

  got <- .gtre_csn_te(e, rep(t_i, N), sig_u, sig_v, sig_r = tiny, sig_h = tiny)

  s2u <- sig_u^2; s2v <- sig_v^2
  eps <- unlist(e)
  want <- .te_battese_coelli(mu_star = -eps * s2u / (s2u + s2v),
                             sigma_star = sig_u * sig_v / sqrt(s2u + s2v))
  expect_equal(length(got$U), length(eps))
  expect_equal(got$U, as.numeric(want), tolerance = 1e-4)
  ## And the test points are not all in a benign interior: one residual is
  ## near zero and one is far into the tail, which is where a quadrature or
  ## orthant-probability error would show up first.
  expect_lt(min(eps), -2.5)
  expect_gt(max(eps), -0.05)
})

test_that("all four scales enter the transient predictor", {
  ## The specific defect: sigma_r and sigma_h were read and then ignored.
  ## Changing sigma_h alone must move U.
  t_i <- 5L; N <- 4L
  set.seed(11L)
  e <- lapply(seq_len(N), function(i) pmin(rnorm(t_i, -0.8, 0.5), 0))

  base <- .gtre_csn_te(e, rep(t_i, N), 0.8, 0.3, sig_r = 0.4, sig_h = 0.8)
  moreh <- .gtre_csn_te(e, rep(t_i, N), 0.8, 0.3, sig_r = 0.4, sig_h = 1.6)
  morer <- .gtre_csn_te(e, rep(t_i, N), 0.8, 0.3, sig_r = 0.9, sig_h = 0.8)

  expect_false(isTRUE(all.equal(base$U, moreh$U)))
  expect_false(isTRUE(all.equal(base$U, morer$U)))
  ## Efficiencies stay in (0, 1].
  for (r in list(base, moreh, morer)) {
    expect_true(all(r$U > 0 & r$U <= 1))
    expect_true(all(r$H > 0 & r$H <= 1))
  }
  expect_equal(length(base$H), N)
  expect_equal(length(base$U), N * t_i)
})

test_that("fiml and sml predict the same transient efficiency on the same data", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 40)
  f_sml  <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = d, model_name = "GTRE",
                                  individual = "name", time = "year",
                                  estimator = "sml"))
  f_fiml <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = d, model_name = "GTRE",
                                  individual = "name", time = "year",
                                  estimator = "fiml"))
  u_s <- as.numeric(f_sml$U); u_f <- as.numeric(f_fiml$U)
  expect_equal(length(u_s), length(u_f))
  ## The two estimators land on slightly different parameter vectors, so their
  ## predictions are not identical -- but they are now the same QUANTITY. The
  ## old code gave a rank correlation of about 0.74 here.
  expect_gt(stats::cor(u_s, u_f), 0.97)
  expect_gt(stats::cor(as.numeric(f_sml$H), as.numeric(f_fiml$H)), 0.97)
})
