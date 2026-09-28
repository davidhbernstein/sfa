## TRE_Z reported a NEGATIVE sigma_r on about 7% of samples.
##
## sigma_r's sign is not identified. It enters the likelihood only as
## `x[2] * R_h1` -- a scale on draws that are symmetric about zero -- so
## +sigma_r and -sigma_r describe the same model and fit equally well. sigma_v
## is different: it also appears in `lamb_fun <- sigma_u_fun / sigma_v_fun`,
## which fixes its sign.
##
## TRE_Z was the only psfm() block that handed bobyqa `lower.bobyqa = -Inf`
## rather than a per-parameter vector flooring the two scales, so it was the
## only one that could reach the negative root -- 110 negative fits across
## 5486 archived TRE_Z replications against 0 for TRE, GTRE, GTRE_Z and
## GTRE_FML. Every later stage floors sigma_r at MIN_POSITIVE, so by the time
## they ran the value was already outside their box and no longer correctable.
##
## The convergence harness read the consequence as non-convergence: sigr's MSE
## was FLAT over n = 1000..5000 (slope +0.010, R2 0.008) while every other
## parameter passed at root-n. On |sigr| the same archive gives slope -0.949.
## The estimator was fine; the sign was not.

test_that("TRE_Z never reports a negative sigma_r", {
  skip_on_cran()
  ## rand = 1 returned sigma_r = -0.2247 before the bound was added, so this
  ## fails loudly if the lower limit is ever relaxed back to -Inf.
  d <- as.data.frame(data_gen_p(t = 10, N = 100, rand = 1, sig_u = 1,
    sig_v = 0.3, sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5,
    beta2 = 0.5))
  f <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
    data = d, individual = "name", keep_objective = TRUE))

  sigr <- unname(coef(f)[2])
  expect_true(is.finite(sigr))
  expect_gte(sigr, 0)
  ## Not merely non-negative: it should be the positive mirror of what the
  ## unbounded run found, near the true 0.2.
  expect_equal(sigr, 0.2238, tolerance = 0.05)

  ## sigma_v's sign IS identified, so it must be positive for a different
  ## reason and by a wide margin.
  expect_gt(unname(coef(f)[1]), 0)

  ## The property that makes the bound necessary rather than cosmetic: the
  ## objective barely moves when sigma_r flips, but moves enormously when
  ## sigma_v does. If this ever stops holding, the rationale above is wrong
  ## and the bound needs rethinking rather than keeping.
  p <- f$opt$par
  base <- f$objective(p)
  flip_r <- f$objective(replace(p, 2L, -p[2L]))
  flip_v <- f$objective(replace(p, 1L, -p[1L]))
  expect_lt(abs(flip_r - base) / abs(base), 1e-3)
  expect_gt(abs(flip_v - base) / abs(base), 1)
})
