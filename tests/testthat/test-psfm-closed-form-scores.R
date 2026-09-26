## Per-firm score contributions for psfm()'s CLOSED-FORM panel likelihoods
## (gap H14 item 1). These sum over firms in closed form rather than by
## simulation, so a contribution is a firm -- not a firm-year -- exactly as in
## the simulated-ML models covered by test-psfm-scores.R.
##
## Two pre-existing defects surfaced while wiring this up and are pinned here:
## PL80's vcov() was returning the estimation-scale matrix under the reported
## names, and PL80_MVTN stored `out` untransposed.

.cf_panel <- function(seed = 21, N = 40, t = 5) {
  as.data.frame(data_gen_p(t = t, N = N, rand = seed, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
}

.cf_check <- function(fit, N, label) {
  par <- fit$opt$par
  ll <- fit$objective(par, per_obs = TRUE)
  expect_length(ll, N)
  expect_true(all(is.finite(ll)), info = label)
  ## The per-firm contributions ARE the objective, un-negated.
  expect_equal(sum(ll), -fit$objective(par), tolerance = 1e-8)
  expect_identical(fit$n_units, N)

  G <- suppressWarnings(estfun.sfareg(fit))
  expect_equal(dim(G), c(N, length(fit$coefficients)), info = label)
  expect_true(all(is.finite(G)), info = label)
  ## First-order conditions: the scores vanish at the optimum. That is the
  ## property that makes them scores rather than arbitrary differences.
  ##
  ## They do NOT hold at a BOUND, where the gradient is free to be non-zero
  ## pointing out of the admissible region, so asserting them there tests a
  ## property the model does not have. estfun() reports whether any difference
  ## left the region; when one did, the assertion is skipped rather than
  ## loosened. FD reaches such a point on the CI platforms and not on macOS,
  ## which is how this surfaced -- see notes/code_history/sandwich_methods.md.
  ## The caller checks that SOME model still exercised the assertion, so this
  ## exemption cannot quietly empty the test.
  foc <- identical(attr(G, "n_bailed"), 0L)
  if (foc) expect_lt(max(abs(colSums(G))) / nrow(G), 1e-3)
  ## bread() scales by the score unit count, not nobs() = N * T.
  expect_equal(suppressWarnings(bread.sfareg(fit)), stats::vcov(fit) * N)
  invisible(foc)
}

test_that("the closed-form time-decay panel models retain per-firm contributions", {
  skip_on_cran()
  d <- .cf_panel()
  asserted <- logical(0)
  for (m in c("PL80", "BC92", "K1990", "K1990modified")) {
    fit <- suppressWarnings(psfm(y_bc92 ~ x1 + x2, model_name = m, data = d,
      individual = "name", time = "year", keep_objective = TRUE))
    asserted <- c(asserted, .cf_check(fit, 40L, m))
  }
  ## Teeth. The boundary exemption in .cf_check() is correct but it must not
  ## be able to skip EVERY model and leave a test that asserts nothing about
  ## the first-order conditions.
  expect_true(any(asserted),
    info = "no closed-form model exercised the FOC assertion")
})

test_that("PL80's vcov() is on the REPORTED scale, not the estimation scale", {
  skip_on_cran()
  d <- .cf_panel()
  fit <- suppressWarnings(psfm(y_bc92 ~ x1 + x2, model_name = "PL80", data = d,
    individual = "name", time = "year", keep_objective = TRUE))

  ## This model ESTIMATES (sigma_v, sigma_u, beta) and REPORTS
  ## (beta, sigmaSq, gamma) -- a permutation as well as a transformation.
  ## Before 1.2.1 vcov() named the estimation-scale inverse Hessian with the
  ## reported names, so confint()'s intercept interval was built from
  ## Var(sigma_v). The Jacobian is stored on the fit precisely to stop that.
  expect_false(is.null(fit$par_scale))
  expect_true(is.matrix(fit$par_scale))
  expect_equal(dim(fit$par_scale),
    c(length(fit$coefficients), length(fit$opt$par)))

  V <- stats::vcov(fit)
  expect_identical(colnames(V), names(fit$coefficients))
  ## The one check that would have caught the bug: vcov() and the printed
  ## standard errors must be the same number.
  expect_equal(unname(sqrt(diag(V))), unname(fit$std.errors), tolerance = 1e-8)

  ## And the raw inverse Hessian must NOT equal it, or the Jacobian is a no-op
  ## and the test above would pass for the wrong reason.
  raw <- sqrt(diag(solve(fit$opt$hessian)))
  expect_false(isTRUE(all.equal(unname(raw), unname(fit$std.errors))))

  ## confint() rides on vcov(), so it inherits the fix.
  ci <- stats::confint(fit)
  expect_equal(unname(ci[, 1]),
    unname(fit$coefficients - stats::qnorm(0.975) * fit$std.errors),
    tolerance = 1e-6)
})

test_that("BHHH and Hessian standard errors agree on the slopes", {
  skip_on_cran()
  d <- .cf_panel()
  fit <- suppressWarnings(psfm(y_bc92 ~ x1 + x2, model_name = "PL80", data = d,
    individual = "name", time = "year", keep_objective = TRUE))
  vb <- sqrt(diag(stats::vcov(fit, type = "bhhh")))
  expect_identical(names(vb), names(fit$coefficients))
  ## Only the slopes. The variance parameters of a time-invariant panel model
  ## are estimated off N firms rather than N*T rows, so the two estimators of
  ## the information matrix diverge much further there.
  sl <- grep("^x[12]$", names(vb), value = TRUE)
  expect_gt(length(sl), 0)
  r <- vb[sl] / fit$std.errors[sl]
  expect_true(all(r > 0.5 & r < 2), info = paste(signif(r, 4), collapse = ", "))
})

test_that("PL80_MVTN stores out as p x 3 like every other entry point", {
  skip_on_cran()
  ## N = 60 rather than 40: at N = 40 this DGP maximises on the sigma_v lower
  ## bound, where the score does not vanish and .cf_check()'s first-order
  ## condition does not apply.
  d <- .cf_panel(seed = 5, N = 60, t = 4)
  fit <- suppressWarnings(psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN",
    data = d, individual = "name", keep_objective = TRUE))
  ## Was 3 x p, so the documented fit$out[, "par"] failed with "subscript out
  ## of bounds" on this model alone.
  expect_equal(dim(fit$out), c(length(fit$coefficients), 3L))
  expect_identical(colnames(fit$out), c("par", "st_err", "t-val"))
  expect_equal(unname(fit$out[, "par"]), unname(fit$coefficients))
  expect_equal(unname(fit$out[, "st_err"]), unname(fit$std.errors))
  .cf_check(fit, 60L, "PL80_MVTN")
})

test_that("PL80_MVTN reaches a maximum rather than stopping short of one", {
  skip_on_cran()
  ## Its likelihood contains an orthant probability computed to limited
  ## precision, so the finite-difference gradient is noisy and a single
  ## L-BFGS-B run stops at points that are not maxima -- while still
  ## reporting convergence = 0. A derivative-free stage was added between two
  ## L-BFGS-B runs in 1.2.1; this pins that the returned point cannot be
  ## improved by perturbing one parameter at a time.
  for (N in c(40L, 60L)) {
    d <- .cf_panel(seed = 5, N = N, t = 4)
    fit <- suppressWarnings(psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN",
      data = d, individual = "name", keep_objective = TRUE))
    par <- fit$opt$par
    v0 <- fit$objective(par)
    lower <- c(1e-3, 1e-3, -1 / (4 - 1) + 0.02, rep(-Inf, length(par) - 3))
    worst <- Inf
    for (j in seq_along(par)) {
      h <- 1e-4 * max(abs(par[j]), 1)
      for (d_ in c(h, -h)) {
        cand <- par
        cand[j] <- cand[j] + d_
        ## A step that leaves the admissible region says nothing about
        ## whether an interior maximum was reached.
        if (cand[j] < lower[j]) next
        worst <- min(worst, fit$objective(cand) - v0)
      }
    }
    expect_gt(worst, -1e-4, label = paste0("N=", N, " directional improvement"))
  }
})

test_that("TFE_WMLE and FD retain per-firm contributions", {
  skip_on_cran()
  d <- .cf_panel()
  fit <- suppressWarnings(psfm(y_tfe ~ x1 + x2, model_name = "TFE_WMLE",
    data = d, individual = "name", time = "year", keep_objective = TRUE))
  .cf_check(fit, 40L, "TFE_WMLE")

  fd <- suppressWarnings(psfm(y_fd ~ x1 + x2 | z_gtre, model_name = "FD",
    data = d, individual = "name", time = "year", keep_objective = TRUE))
  .cf_check(fd, 40L, "FD")

  ## FD's Hessian is routinely indefinite, so several of its printed standard
  ## errors come back NaN. The OPG needs no Hessian, so it is defined exactly
  ## where that path fails -- which is the point of offering it.
  vb <- sqrt(diag(stats::vcov(fd, type = "bhhh")))
  expect_true(all(is.finite(vb)))
})

test_that("GTRE_FML retains per-firm contributions", {
  skip_on_cran()
  ## A sample with a genuine persistent component. When sigr collapses to the
  ## zero boundary the numerical score is meaningless, because perturbing it
  ## downwards leaves the admissible region and hits the likelihood penalty --
  ## that is a property of the boundary, not of the score code.
  d <- as.data.frame(data_gen_p(t = 5, N = 40, rand = 3, sig_u = 0.8,
    sig_v = 0.4, sig_r = 0.7, sig_h = 0.6, cons = 1, beta1 = 0.5, beta2 = 0.5))
  fit <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE_FML",
    data = d, individual = "name", keep_objective = TRUE))
  skip_if(isTRUE(fit$sigma_r_at_bound), "sigr collapsed in this sample")
  .cf_check(fit, 40L, "GTRE_FML")
})

test_that("clustered standard errors now work for the closed-form panel models", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  d <- .cf_panel()
  fit <- suppressWarnings(psfm(y_bc92 ~ x1 + x2, model_name = "PL80", data = d,
    individual = "name", time = "year", keep_objective = TRUE))
  ## Contributions are per firm, so the cluster variable is one entry per
  ## FIRM, in firm order -- a coarser grouping such as a region.
  set.seed(4)
  region <- sample(seq_len(8), fit$n_units, replace = TRUE)
  V <- sandwich::vcovCL(fit, cluster = region)
  expect_equal(dim(V), c(length(fit$coefficients), length(fit$coefficients)))
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))
  expect_true(all(is.finite(sandwich::sandwich(fit))))

  ## A firm-year-length cluster vector is the natural mistake, and sandwich
  ## catches it rather than silently recycling.
  expect_error(sandwich::vcovCL(fit, cluster = d$name), "cluster")
})
