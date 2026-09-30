## Issue #30: from its single moment start (shape 0.5), NNAK's stages could
## drift up in shape and stop far below the maximum -- below even the
## half-normal NNAK nests at shape 0.5. On the sample below main returned
## logLik -428.92 at shape 2.34; NHN gives -404.89; the profile over the shape
## (computed independently in the issue) peaks at -359.55 near shape 0.005.

nnak_d <- function(seed) {
  as.data.frame(data_gen_cs(N = 200, rand = seed, cons = 0.5, beta1 = 0.5,
    beta2 = 0.5, sig_u = 1, sig_v = 1, mu = 0.5, a = 5))
}

test_that("NNAK at shape 0.5 is exactly the half-normal it is compared to", {
  skip_on_cran()
  ## The bound the fit is held to below rests on this identity, so check it
  ## rather than assume it. At shape 1/2, u^2 ~ Gamma(1/2, scale 2 sigma_u^2),
  ## i.e. sigma_u^2 times a chi-square(1), so u = sigma_u |Z|: the half-normal.
  d <- nnak_d(1002412)
  k <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NNAK",
    keep_objective = TRUE))
  h <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NHN"))
  hp <- h$out[, "par"]
  lam <- hp[["lambda"]]; sig <- hp[["sigma"]]
  q <- c(sig / sqrt(1 + lam^2), lam * sig / sqrt(1 + lam^2), 0.5, hp[-(1:2)])
  expect_equal(-k$objective(q), as.numeric(logLik(h)), tolerance = 1e-9)
})

test_that("NNAK no longer stops below the half-normal it nests (issue #30)", {
  skip_on_cran()
  d <- nnak_d(1002412)
  k <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NNAK"))
  h <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NHN"))
  expect_gte(as.numeric(logLik(k)), as.numeric(logLik(h)) - 1e-6)
  ## And it reaches the profile maximum rather than merely clearing NHN.
  expect_gt(as.numeric(logLik(k)), -359.6)
  expect_lt(k$out["mu", "par"], 0.05)
  expect_true(all(is.finite(k$exp_u_hat)))
})

test_that("NNAK reports its starting-value search", {
  skip_on_cran()
  k <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = nnak_d(1002412),
    model_name = "NNAK"))
  expect_false(is.null(k$ng_starts))
  expect_gte(k$ng_starts$n_tried, 5L)
  ## The fit is never worse than the best polished start (gap A25's guard).
  expect_gte(-k$opt$value, k$ng_starts$best - 1e-6)
})

test_that("the NNAK candidates put E[u] where the moments do", {
  ## Each candidate's E[u] = sigma_u Gamma(m + 1/2) / (Gamma(m) sqrt(m)) either
  ## equals the moment anchor or, where that would over-spend the variance,
  ## leaves sigma_v a positive share. Checked by simulation, not by the
  ## formula the helper uses.
  set.seed(3)
  e <- rnorm(2000, 0, 0.5) - abs(rnorm(2000, 0, 1))
  cands <- .nnak_start_candidates(e, beta_0_st = 1, beta_hat = 0.5)
  expect_gte(length(cands), 5L)
  for (z in cands) {
    sv <- z[1]; su <- z[2]; m <- z[3]
    expect_gt(sv, 0)
    u <- sqrt(stats::rgamma(2e5, shape = m, scale = su^2 / m))
    expect_equal(z[4], 1 + mean(u), tolerance = 0.02)
    expect_lte(stats::var(u) + sv^2, stats::var(e) * 1.05)
  }
})
