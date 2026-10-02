## Issue #55, defect 1: at the lambda floor the NTN log-likelihood was formed
## as log Phi(aa) - log Phi(bb) with a shared mu / lam inside both arguments.
## With mu < 0 each term reaches O(-1e13) while the difference is O(1), so the
## answer was rounding noise -- it came out ABOVE the true value and above the
## OLS maximum. See notes/code_history/sfm.md.

gen <- function(seed) {
  as.data.frame(data_gen_cs(
    N = 200, rand = seed, cons = 0.5, beta1 = 0.5, beta2 = 0.5,
    sig_u = 1, sig_v = 1, mu = 0.5, a = 5
  ))
}

## As lam -> 0 the NTN density collapses to a normal in eps: l3 + l4 + l5 ->
## -eps^2 (1 + lam^2) / (2 sigma^2) because aa / bb -> 1. For mu < 0 the mean
## is 0, for mu >= 0 it is -mu. Independent of the package's own algebra.
lam0_limit <- function(y, X, beta, sigma, mu) {
  eps <- as.numeric(y - X %*% beta)
  sum(stats::dnorm(eps + if (mu < 0) 0 else mu, 0, sigma, log = TRUE))
}

test_that("NTN reaches its lambda -> 0 limit instead of cancelling (#55)", {
  d <- gen(1018875)
  X <- cbind(1, d$x1, d$x2)
  ols <- stats::lm(y_pcs_ge ~ x1 + x2, data = d)
  f <- suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", keep_objective = TRUE
  ))
  expect_true(is.function(f$objective))

  ## A fixed parameter vector, so this tests the likelihood and not the
  ## optimizer: OLS betas and sigma, mu on the negative side where the old
  ## form cancelled.
  sig <- stats::sd(stats::resid(ols))
  mk <- function(lam, mu) c(lam, sig, mu, unname(stats::coef(ols)))

  for (mu in c(-1.0061, -0.5, -2)) {
    target <- lam0_limit(d$y_pcs_ge, X, stats::coef(ols), sig, mu)
    ## The gap to the limit is O(lam^2), so it must SHRINK as lam falls. On the
    ## old form it turned around below lam = 1e-4 and grew without bound.
    errs <- vapply(c(1e-4, 1e-5, 1e-6, 1e-7), function(lam) {
      abs(-f$objective(mk(lam, mu)) - target)
    }, numeric(1))
    expect_true(all(diff(errs) < 0))
    expect_lt(errs[length(errs)], 1e-8)
  }
})

test_that("NTN cannot beat OLS at the lambda floor (#55)", {
  ## The lam -> 0 limit is a Gaussian log-likelihood at a (beta, sigma) that
  ## OLS maximises over, so near the floor NTN is bounded above by OLS. The old
  ## form reported 1.21 ABOVE it on this sample.
  for (cs in list(
    list(seed = 1018875, col = "y_pcs_ge"),
    list(seed = 1013522, col = "y_pcs_thn"),
    list(seed = 1011401, col = "y_pcs_st")
  )) {
    d <- gen(cs$seed)
    fo <- stats::as.formula(paste(cs$col, "~ x1 + x2"))
    ll_ols <- as.numeric(stats::logLik(stats::lm(fo, data = d)))
    f <- suppressWarnings(sfm(fo, data = d, model_name = "NTN", keep_objective = TRUE))
    ols <- stats::lm(fo, data = d)
    sig <- stats::sd(stats::resid(ols))
    for (lam in c(1e-5, 1e-6, 1e-7)) {
      for (mu in c(-0.25, -0.8544, -1.5)) {
        p <- c(lam, sig, mu, unname(stats::coef(ols)))
        expect_lte(-f$objective(p), ll_ols + 1e-6)
      }
    }
    ## And the fit the optimizer actually returns obeys it too when it lands
    ## on the floor.
    p_hat <- f$out[, "par"]
    if (p_hat[["lambda"]] < 1e-5) {
      expect_lte(as.numeric(stats::logLik(f)), ll_ols + 1e-6)
    }
  }
})

test_that("the stable NTN form is unchanged away from the floor (#55)", {
  ## Where log Phi(aa) - log Phi(bb) does not cancel, the rewrite must be the
  ## same number: this reference is the pre-#55 expression, written out here.
  d <- gen(1018875)
  X <- cbind(1, d$x1, d$x2)
  ols <- stats::lm(y_pcs_ge ~ x1 + x2, data = d)
  f <- suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", keep_objective = TRUE
  ))
  beta <- unname(stats::coef(ols))
  sig <- stats::sd(stats::resid(ols))
  eps <- as.numeric(d$y_pcs_ge - X %*% beta)

  reference <- function(lam, mu) {
    sum(-log(sig^2) / 2 - log(2 * pi) / 2 -
      (1 / (2 * sig^2)) * (-eps - mu)^2 +
      stats::pnorm(((mu / lam) - eps * lam) / sig, log.p = TRUE) -
      stats::pnorm((mu / sig) * sqrt(1 + lam^(-2)), log.p = TRUE))
  }
  for (lam in c(2, 1, 0.5, 0.1, 0.01, 1e-3)) {
    for (mu in c(-1.5, -0.5, 0.5, 1.5)) {
      expect_equal(-f$objective(c(lam, sig, mu, beta)), reference(lam, mu),
        tolerance = 1e-7
      )
    }
  }
})

test_that("NTN with mu > 0 at the floor still takes the direct path (#55)", {
  ## aa and bb both go to +Inf there, log Phi of each goes to 0, and nothing
  ## cancels -- so that branch must keep using the unmodified expression.
  d <- gen(1018875)
  X <- cbind(1, d$x1, d$x2)
  ols <- stats::lm(y_pcs_ge ~ x1 + x2, data = d)
  f <- suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", keep_objective = TRUE
  ))
  beta <- unname(stats::coef(ols))
  sig <- stats::sd(stats::resid(ols))
  eps <- as.numeric(d$y_pcs_ge - X %*% beta)
  for (lam in c(1e-3, 1e-5, 1e-7)) {
    for (mu in c(0.5, 2)) {
      direct <- sum(-log(sig^2) / 2 - log(2 * pi) / 2 -
        (1 / (2 * sig^2)) * (-eps - mu)^2 +
        stats::pnorm(((mu / lam) - eps * lam) / sig, log.p = TRUE) -
        stats::pnorm((mu / sig) * sqrt(1 + lam^(-2)), log.p = TRUE))
      expect_equal(-f$objective(c(lam, sig, mu, beta)), direct, tolerance = 1e-12)
    }
  }
})

## Issue #55, defect 2: NTN at mu = 0 IS NHN with the same (lambda, sigma), and
## mu = 0 is interior -- start_cs() bounds only lambda and sigma below, leaving
## mu free -- so NTN's maximum cannot lie below NHN's. It did, by up to 15.9.
## sfm() now checks the finished fit against a mu = 0 fit and polishes from it,
## the way NG checks against shape 1 (#45) and NNAK against m = 0.5 (#30).

test_that("mu = 0 is interior for NTN, which is what makes NHN a floor (#55)", {
  ## The likelihood is finite and smooth on BOTH sides of mu = 0, so mu = 0 is
  ## an interior point and not a bound. That is what separates this comparison
  ## from a boundary supremum such as the sigma_v collapse or TSL's OLS limit:
  ## a maximum below NHN's is a failure, not a limit that is merely approached.
  d <- gen(1018875)
  ols <- stats::lm(y_pcs_ge ~ x1 + x2, data = d)
  f <- suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", keep_objective = TRUE
  ))
  sig <- stats::sd(stats::resid(ols))
  at <- function(mu) -f$objective(c(1.5, sig, mu, unname(stats::coef(ols))))
  vals <- vapply(c(-0.5, -0.05, -1e-5, 0, 1e-5, 0.05, 0.5), at, numeric(1))
  expect_true(all(is.finite(vals)))
  ## Continuous through zero from both sides, to the step size.
  expect_lt(abs(at(-1e-6) - at(0)), 1e-3)
  expect_lt(abs(at(1e-6) - at(0)), 1e-3)
  ## And a fit is free to return mu on either side of zero.
  expect_gt(at(-0.5), -Inf)
})

test_that("NTN does not end below the nested NHN (#55)", {
  for (cs in list(
    list(seed = 1018875, col = "y_pcs_ge"),
    list(seed = 1005038, col = "y_pcs_tsl"),
    list(seed = 1011401, col = "y_pcs_st"),
    list(seed = 1013522, col = "y_pcs_thn")
  )) {
    d <- gen(cs$seed)
    fo <- stats::as.formula(paste(cs$col, "~ x1 + x2"))
    ll_ntn <- as.numeric(stats::logLik(suppressWarnings(
      sfm(fo, data = d, model_name = "NTN")
    )))
    ll_nhn <- as.numeric(stats::logLik(suppressWarnings(
      sfm(fo, data = d, model_name = "NHN")
    )))
    expect_gte(ll_ntn, ll_nhn - 1e-6)
  }
})

test_that("the mu = 0 reference is a real NHN fit, not just a start (#55)", {
  ## Held at mu = 0 the NTN likelihood must agree with NHN's own at the same
  ## (lambda, sigma, beta) -- that identity is what licenses the comparison.
  d <- gen(1018875)
  fo <- y_pcs_ge ~ x1 + x2
  nhn <- suppressWarnings(sfm(fo, data = d, model_name = "NHN"))
  f <- suppressWarnings(sfm(fo, data = d, model_name = "NTN", keep_objective = TRUE))
  p_nhn <- nhn$out[, "par"]
  lifted <- c(p_nhn[["lambda"]], p_nhn[["sigma"]], 0, unname(p_nhn[3:5]))
  expect_equal(-f$objective(lifted), as.numeric(stats::logLik(nhn)),
    tolerance = 1e-8
  )
})

test_that("the mu = 0 reference copes with every NTN parameter layout (#55)", {
  ## .nhn_ref is built by inserting mu = 0 into a vector whose length varies
  ## with uhet/muhet and with the intercept, and one of its seeds comes from
  ## start_v_nhn, which carries no het coefficients at all. Each of these fits
  ## must still complete and report the right number of parameters.
  d <- gen(1018875)
  d$z1 <- d$z
  d$z2 <- d$zp
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN"))$out), 6L)
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", uhet = ~z1))$out), 7L)
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", muhet = ~z1))$out), 7L)
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ x1 + x2,
    data = d, model_name = "NTN", uhet = ~z1, muhet = ~z2))$out), 8L)
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ 0 + x1 + x2,
    data = d, model_name = "NTN"))$out), 5L)
  ## A numeric start_val skips the reference entirely, as it does for NG.
  expect_equal(nrow(suppressWarnings(sfm(y_pcs_ge ~ x1 + x2, data = d,
    model_name = "NTN", start_val = c(1.5, 2, 0, 0.5, 0.5, 0.5)))$out), 6L)
})
