## Issue #28: TSL's composed log-density formed each tilt,
## sig_v^2 / (2 sig_u^2) + eps / sig_u, separately from log Phi(aa). The two
## are each about sig_v^2 / (2 sig_u^2) with opposite signs, so once sigma_u
## approaches its floor the sum is rounding noise -- and the optimizer climbed
## the positive noise to a reported logLik of +1.7e40.

test_that("TSL's likelihood tends to its normal limit as sigma_u collapses", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  f <- suppressWarnings(sfm(y_pcs_tsl ~ x1 + x2, model_name = "TSL", data = d,
    keep_objective = TRUE))
  p <- f$out[, "par"]
  X <- cbind(1, d$x1, d$x2)
  for (su in c(1e-7, 1e-9, 1e-12, 1e-16)) {
    q <- p
    q[2] <- su
    pc <- as.numeric(f$objective(q, per_obs = TRUE))
    ## u degenerates to 0 as sigma_u -> 0, leaving eps ~ N(0, sigma_v^2). The
    ## first omitted term is O(sigma_u / sigma_v), below 1e-5 here.
    lim <- dnorm(as.numeric(d$y_pcs_tsl - X %*% q[4:6]), 0, q[1], log = TRUE)
    expect_equal(pc, lim, tolerance = 1e-5, info = su)
  }
})

test_that("the point the optimizer used to converge to scores as what it is", {
  skip_on_cran()
  ## The reproduction in issue #28: sfm() returned sigma_v ~ 1.13e20 with
  ## sigma_u on its floor and reported logLik 1.7e40. At that point u is
  ## negligible against sigma_v, so the density is the normal one and the
  ## log-likelihood is hugely negative, not positive.
  d <- as.data.frame(data_gen_cs(N = 200, rand = 1002715, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, sig_u = 1, sig_v = 1, mu = 0.5, a = 5))
  f <- suppressWarnings(sfm(y_pcs_thn ~ x1 + x2, model_name = "TSL", data = d,
    keep_objective = TRUE))
  q <- f$out[, "par"]
  q[1:3] <- c(1.13e20, 1e-7, 6.53e6)
  pc <- as.numeric(f$objective(q, per_obs = TRUE))
  X <- cbind(1, d$x1, d$x2)
  lim <- dnorm(as.numeric(d$y_pcs_thn - X %*% q[4:6]), 0, q[1], log = TRUE)
  expect_true(all(pc < 0))
  expect_equal(pc, lim, tolerance = 1e-8)
  ## And the fit itself is no worse than OLS, TSL's sigma_u -> 0 limit.
  expect_gte(as.numeric(logLik(f)),
    as.numeric(logLik(stats::lm(y_pcs_thn ~ x1 + x2, data = d))) - 1e-6)
})

test_that("the tilt form equals the old form where the old form is exact", {
  ## At ordinary parameter values the rewrite must change nothing.
  old <- function(e, sv, su, lam) {
    A <- sv^2 / (2 * su^2) + e / su
    B <- (1 + lam)^2 * sv^2 / (2 * su^2) + e * (1 + lam) / su
    l1 <- log(2) + A + pnorm(-sv / su - e / sv, log.p = TRUE)
    l2 <- B + pnorm(-sv * (1 + lam) / su - e / sv, log.p = TRUE)
    log1p(lam) - log(2 * lam + 1) - log(su) +
      l1 + log(-expm1(pmin(l2 - l1, -.Machine$double.eps)))
  }
  new <- function(e, sv, su, lam) {
    q <- -e^2 / (2 * sv^2)
    l1 <- log(2) + q + .log_phi_tilt(-sv / su - e / sv)
    l2 <- q + .log_phi_tilt(-sv * (1 + lam) / su - e / sv)
    log1p(lam) - log(2 * lam + 1) - log(su) +
      l1 + log(-expm1(pmin(l2 - l1, -.Machine$double.eps)))
  }
  g <- expand.grid(e = c(-3, -1, 0, 0.5, 2), sv = c(0.3, 1), su = c(0.2, 0.8),
    lam = c(0.5, 2))
  expect_equal(mapply(new, g$e, g$sv, g$su, g$lam),
    mapply(old, g$e, g$sv, g$su, g$lam), tolerance = 1e-12)
})
