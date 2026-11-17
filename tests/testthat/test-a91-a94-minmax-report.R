## Gaps A91-A94, reported 2026-11-17 from an outside study that fitted these
## models at scale against an independent quadrature likelihood. Each test pins
## arithmetic or an invariant, never where an optimizer lands.

test_that(".log_exp_phi() is exact where the old tilt form cancelled (A92)", {
  ## e < 0 with sigma_v -> 0: -e^2/(2 sv^2) and z^2/2 are each ~7e12 and the
  ## old form lost ~1e-3 here. The closed form is c e + c^2 sv^2/2 + log Phi(z),
  ## and log Phi(z) is 0 to double precision at z ~ 4e6.
  e <- c(-1.423087, -0.835742); sv <- 3.844434e-07; cc <- 1 / 0.1757172
  expect_equal(.log_exp_phi(e, sv, cc), cc * e + cc^2 * sv^2 / 2, tolerance = 1e-14)
  ## Continuous across the switch: both branches agree on either side of it.
  sv2 <- 0.3; cc2 <- 2
  for (z in .SFA_CONSTANTS$NE_TILT_SWITCH + c(-1e-6, 1e-6)) {
    e2 <- -(z + cc2 * sv2) * sv2
    lo <- -e2^2 / (2 * sv2^2) + .log_phi_tilt(z)
    hi <- cc2 * e2 + cc2^2 * sv2^2 / 2 + stats::pnorm(z, log.p = TRUE)
    expect_equal(.log_exp_phi(e2, sv2, cc2), lo, tolerance = 1e-12)
    expect_equal(lo, hi, tolerance = 1e-12)
  }
})

test_that("sfm('NGE') log-likelihood is exact near sigma_v -> 0 (A92)", {
  skip_on_cran()
  set.seed(336)
  n <- 60; x1 <- rnorm(n); x2 <- rnorm(n)
  lam <- 1 / 0.18
  u <- -log(1 - sqrt(runif(n))) / lam   # GE(2, lambda): max of two exponentials
  d <- data.frame(y = 0.5 + 0.5 * x1 + 0.5 * x2 - u, x1 = x1, x2 = x2)
  f <- suppressWarnings(sfm(y ~ x1 + x2, data = d, model_name = "NGE",
                            keep_objective = TRUE))
  ## Evaluated at a FIXED point with sigma_v = 4e-7, not at the fit.
  p <- c(4e-7, 0.18, 0.5, 0.5, 0.5)
  e <- d$y - cbind(1, x1, x2) %*% p[3:5]
  ref <- sum(vapply(as.numeric(e), function(ei) {
    g <- function(u) dnorm(ei + u, 0, p[1]) * 2 * lam * exp(-lam * u) * (-expm1(-lam * u))
    at <- max(0, -ei)
    br <- sort(unique(c(0, if (at > 0) c(max(0, at - 50 * p[1]), at, at + 50 * p[1]), 50 * p[1], Inf)))
    log(sum(vapply(seq_len(length(br) - 1L), function(k)
      integrate(g, br[k], br[k + 1L], rel.tol = 1e-12, subdivisions = 5000L)$value, 0)))
  }, 0))
  expect_equal(-f$objective(p), ref, tolerance = 1e-8)
})

test_that("cw_test() refuses NU, whose moments cannot identify it (A93)", {
  set.seed(1)
  e <- rnorm(300, 0, 0.3) - runif(300, 0, 1)
  expect_error(cw_test(e, model_name = "NU"), "uniform inefficiency is symmetric")
})

test_that("cw_test('NHN') stops cleanly past the half-normal skewness limit (A94)", {
  set.seed(2)
  e <- rnorm(400, 0, 0.05) - rlnorm(400, 0, 1)   # skewness well past 0.995
  ec <- e - mean(e)
  expect_lt(mean(ec^3) / mean(ec^2)^1.5, -0.995)
  expect_error(cw_test(e, model_name = "NHN"), "cannot reproduce")
  ## Half-normal data still gives a statistic.
  e2 <- rnorm(400, 0, 0.3) - abs(rnorm(400, 0, 1))
  o <- cw_test(e2, model_name = "NHN")
  expect_true(is.finite(o$statistic) && o$p.value > 0 && o$p.value <= 1)
})

test_that("sfm('NW') starts include near-noise-free points and k is capped (A91)", {
  e <- rnorm(200, 0, 0.3) - rweibull(200, 3, 0.7)
  cands <- .nw_start_candidates(e, 0, c(0.5, 0.5))
  sd_e <- sqrt(mean((e - mean(e))^2))
  ks <- vapply(cands, `[`, numeric(1), 3)
  svs <- vapply(cands, `[`, numeric(1), 1)
  expect_true(all(c(1, 1.5, 2, 3, 5) %in% ks))
  expect_true(any(abs(svs - 0.05 * sd_e) < 1e-12))
  expect_true(all(ks <= .NW_K_RIDGE))
})
