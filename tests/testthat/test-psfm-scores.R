## Per-firm score contributions for psfm()'s simulated-ML panel models (gap A5).
## Before 1.2.1 these closures had no `per_obs` branch, so keep_objective = TRUE
## stored an objective that estfun(), sandwich(), vcov(type = "bhhh") and
## influence_sfa() could not use.

.score_panel <- function(seed = 21) {
  as.data.frame(data_gen_p(t = 5, N = 40, rand = seed, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
}

.check_panel_scores <- function(fit, N) {
  par <- fit$opt$par
  ll <- fit$objective(par, per_obs = TRUE)
  expect_length(ll, N)
  expect_true(all(is.finite(ll)))
  ## The per-firm contributions ARE the objective, un-negated.
  expect_equal(sum(ll), -fit$objective(par), tolerance = 1e-10)
  expect_identical(fit$n_units, N)
  G <- estfun.sfareg(fit)
  expect_equal(dim(G), c(N, length(fit$coefficients)))
  ## bread() must scale by the score unit count, not nobs() = N * T.
  expect_equal(bread.sfareg(fit), stats::vcov(fit) * N)
  V <- stats::vcov(fit, type = "bhhh")
  expect_true(all(is.finite(V)))
  if (requireNamespace("sandwich", quietly = TRUE)) {
    expect_true(all(is.finite(sandwich::sandwich(fit))))
  }
  expect_equal(nrow(influence_sfa(fit)$influence), N)
}

test_that("TRE_Z retains per-firm log-likelihood contributions", {
  skip_on_cran()
  d <- .score_panel()
  fit <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
    data = d, individual = "name", time = "year", halton_num = 50,
    keep_objective = TRUE))
  .check_panel_scores(fit, 40L)
})

test_that("GTRE_Z retains per-firm log-likelihood contributions", {
  skip_on_cran()
  d <- .score_panel()
  fit <- suppressWarnings(psfm(y_gtre_zz ~ x1 + x2 | z_gtre | zp_gtre, model_name = "GTRE_Z",
    data = d, individual = "name", time = "year", halton_num = 50,
    keep_objective = TRUE))
  .check_panel_scores(fit, 40L)
})

test_that("GTRE (sml): per-firm contributions, and the loop and vectorised paths agree", {
  skip_on_cran()
  d <- .score_panel()
  fit <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE", estimator = "sml",
    data = d, individual = "name", time = "year", halton_num = 50,
    keep_objective = TRUE))
  .check_panel_scores(fit, 40L)
  par <- fit$opt$par
  old <- options(sfa.gtre_vectorized = FALSE)
  on.exit(options(old), add = TRUE)
  loop <- fit$objective(par, per_obs = TRUE)
  options(sfa.gtre_vectorized = TRUE)
  vec <- fit$objective(par, per_obs = TRUE)
  expect_equal(loop, vec, tolerance = 1e-8)
})

test_that("a psfm() fit without a per-observation likelihood still says so", {
  skip_on_cran()
  d <- .score_panel()
  fit <- suppressWarnings(psfm(y_fd ~ x_fd | z_fd, model_name = "FD", data = d,
    individual = "name", time = "year"))
  expect_error(estfun.sfareg(fit), "does not retain its likelihood")
})
