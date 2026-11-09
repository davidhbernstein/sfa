## A27 (NR) and A59 (NTN): the likelihood's supremum can be the sigma_v -> 0
## corner while the optimizer stops at an interior point below it. The fit now
## probes the corner itself; these tests hold it to that, whichever side of the
## corner each sample's fit ends on.

test_that("an NR fit is not beaten by its own corner probe, and a corner fit says so", {
  skip_on_cran()
  for (s in c(1, 3, 8)) {
    d <- as.data.frame(data_gen_cs(N = 500, rand = s, sig_u = 1, sig_v = 0.1,
      cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5))
    f <- suppressWarnings(sfm(y_pcs_r ~ x1 + x2, model_name = "NR", data = d,
      keep_objective = TRUE))
    p <- .nr_corner_probe(f$objective, f$opt$par, cbind(1, d$x1, d$x2), d$y_pcs_r, 1)
    if (!is.null(p)) expect_gte(as.numeric(logLik(f)), -p$value - 1e-4)
    expect_type(f$sigma_v_at_bound, "logical")
    if (isTRUE(f$sigma_v_at_bound)) expect_true(all(is.na(f$out[, "st_err"])))
  }
})

test_that("an NTN fit is not beaten by its own lambda ladder", {
  skip_on_cran()
  ## seed 3 at sig_v = 0.5 is A59's case: interior lambda 3.5, 1.05 short.
  for (s in c(1, 3, 10)) {
    d <- as.data.frame(data_gen_cs(N = 150, rand = s, cons = 0.5, beta1 = 0.5,
      beta2 = 0.5, sig_u = 1, sig_v = 0.5, mu = 0.5, a = 5))
    f <- suppressWarnings(sfm(y_pcs_tn ~ x1 + x2, model_name = "NTN", data = d,
      keep_objective = TRUE))
    l <- .ntn_lambda_ladder(f$objective, f$opt$par)
    if (!is.null(l)) expect_gte(as.numeric(logLik(f)), -l$value - 1e-4)
    if (isTRUE(f$ntn_sigma_v_at_bound)) expect_true(all(is.na(f$out[, "st_err"])))
  }
})

test_that("the NTN corner the ladder reports is the exact sigma_v = 0 likelihood", {
  skip_on_cran()
  d <- as.data.frame(data_gen_cs(N = 150, rand = 3, cons = 0.5, beta1 = 0.5,
    beta2 = 0.5, sig_u = 1, sig_v = 0.5, mu = 0.5, a = 5))
  f <- suppressWarnings(sfm(y_pcs_tn ~ x1 + x2, model_name = "NTN", data = d))
  skip_if_not(isTRUE(f$ntn_sigma_v_at_bound))
  p <- f$out[, "par"]
  e <- d$y_pcs_tn - cbind(1, d$x1, d$x2) %*% p[4:6]
  expect_true(all(e <= 1e-6))
  ## truncated normal u = -e, the lambda -> Inf limit of the NTN density
  lim <- sum(dnorm((-e - p[["mu"]]) / p[["sigma"]], log = TRUE) - log(p[["sigma"]])) -
    length(e) * pnorm(p[["mu"]] / p[["sigma"]], log.p = TRUE)
  expect_equal(as.numeric(logLik(f)), lim, tolerance = 1e-6)
})
