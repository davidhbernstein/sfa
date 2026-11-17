## Gap A90: copsfm() started sigma_v from Var(eps) = Var(v) + Var(u), floored at
## 0.001. Under dependence that difference goes negative, and the floor put the
## optimizer in the basin of a local maximum at sigma_v -> 0: 11 of 25 fits on
## copula data ended 11 to 129 log-units below the true parameters' likelihood.
## The start rule is tested directly -- it is what the fix guarantees -- rather
## than where any fit lands.

test_that("the sigma_v start is the moment estimate where the moments identify it (A90)", {
  ## m2 = 0.09 + (1 - 2/pi) * 1: exactly sigma_v = 0.3 beside sigma_u = 1.
  m2 <- 0.3^2 + (1 - 2 / pi) * 1^2
  expect_equal(.cop_sv_start(m2, 1, sqrt(m2), "hnormal", "normal"), 0.3)
  ## Laplace noise: Var(v) = 2 delta^2.
  m2l <- 2 * 0.3^2 + (1 - 2 / pi)
  expect_equal(.cop_sv_start(m2l, 1, sqrt(m2l), "hnormal", "laplace"), 0.3)
})

test_that("a degenerate moment estimate does not start sigma_v on a floor (A90)", {
  ## Dependence makes m2 - c * su0^2 negative.
  sd_e <- 0.6
  s0 <- .cop_sv_start(0.2, 1, sd_e, "hnormal", "normal")
  expect_equal(s0, 0.5 * sd_e)
  expect_gt(s0, 0.1)
  ## Generalized logistic: scaled by its variance constant pi^2/3.
  expect_equal(.cop_sv_start(0.2, 1, sd_e, "hnormal", "glogistic"), 0.5 * sd_e / sqrt(pi^2 / 3))
  ## A positive but negligible difference counts as unidentified too.
  tiny <- (1 - 2 / pi) + 1e-6
  expect_equal(.cop_sv_start(tiny, 1, sd_e, "hnormal", "normal"), 0.5 * sd_e)
})
