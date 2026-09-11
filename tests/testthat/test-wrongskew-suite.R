## The wrong-skewness suite: ACOLS (Parmeter and Zhao 2023), constrained MLE
## (Zhao and Parmeter 2022) and the third-moment decomposition of Bonanno,
## De Giovanni and Domma (2017).

## --- ACOLS ------------------------------------------------------------------

test_that("the absolute-moment equations are the ones in the paper", {
  ## Eq. (8): (1/n) sum |sqrt(pi/2) e_i - sigma_u| - (r2 + (2/pi) sigma_u^2)^(1/2)
  set.seed(1)
  e <- rnorm(200); e <- e - mean(e); r2 <- mean(e^2)
  for (su in c(0, 0.3, 1.2)) {
    expect_equal(sfa:::.abs_mom_nhn(su, e, r2),
      mean(abs(sqrt(pi / 2) * e - su)) - sqrt(r2 + (2 / pi) * su^2))
  }
  ## Eq. (13). The 2 exp(lam^2/2) Phi(-lam) term is written through
  ## log Phi so that large lambda does not overflow, which is the reason the
  ## paper had to stop at sigma_u = 0.03.
  for (lam in c(0.1, 1, 5)) {
    expect_equal(sfa:::.abs_mom_ne(lam, e, r2),
      mean(abs(e * sqrt((1 + lam^2) / r2) - 1)) -
        2 * exp(lam^2 / 2) * pnorm(-lam) - sqrt(2 / pi) * lam,
      tolerance = 1e-10, info = lam)
  }
  ## ... and at lambda = 60 the naive form IS Inf while ours is finite.
  expect_false(is.finite(2 * exp(60^2 / 2) * pnorm(-60)))
  expect_true(is.finite(sfa:::.abs_mom_ne(60, e, r2)))
})

test_that("ACOLS solves its own moment conditions", {
  set.seed(2)
  n <- 300
  x <- rnorm(n)
  y <- 1 + 0.8 * x + rnorm(n, 0, 0.6) - abs(rnorm(n, 0, 1))
  X <- cbind("(Intercept)" = 1, x = x)
  f <- sfa:::.cols_fit(y, X, "NHN", 1L, moment = "absolute")
  e <- lm.fit(X, y)$residuals; e <- e - mean(e); r2 <- mean(e^2)
  expect_lt(abs(sfa:::.abs_mom_nhn(f$sigma_u, e, r2)), 1e-8)
  ## and the second condition holds by construction
  expect_equal(f$sigma_v^2 + f$sigma_u^2 * (pi - 2) / pi, r2, tolerance = 1e-8)
  ## the third moment is NOT used: the fit is unchanged if we only knew m3
  expect_identical(f$moment, "absolute")

  y2 <- 1 + 0.8 * x + rnorm(n, 0, 0.6) - rexp(n, 1)
  f2 <- sfa:::.cols_fit(y2, X, "NE", 1L, moment = "absolute")
  e2 <- lm.fit(X, y2)$residuals; e2 <- e2 - mean(e2); r22 <- mean(e2^2)
  expect_lt(abs(sfa:::.abs_mom_ne(f2$sigma_v / f2$sigma_u, e2, r22)), 1e-8)
  expect_equal(f2$sigma_v^2 + f2$sigma_u^2, r22, tolerance = 1e-8)
})

test_that("ACOLS keeps the OLS slopes and corrects only the intercept", {
  set.seed(3)
  n <- 200
  x1 <- rnorm(n); x2 <- runif(n)
  y <- 2 + x1 - 0.5 * x2 + rnorm(n, 0, 0.5) - abs(rnorm(n, 0, 0.9))
  X <- cbind("(Intercept)" = 1, x1 = x1, x2 = x2)
  ols <- lm.fit(X, y)$coefficients
  a <- sfa:::.cols_fit(y, X, "NHN", 1L, moment = "absolute")
  c3 <- sfa:::.cols_fit(y, X, "NHN", 1L, moment = "third")
  expect_equal(unname(a$beta[-1]), unname(ols[-1]))
  expect_equal(unname(c3$beta[-1]), unname(ols[-1]))
  ## Same slopes, DIFFERENT intercept: the two estimators differ only through
  ## the E[u] shift they infer.
  expect_false(isTRUE(all.equal(a$beta[[1]], c3$beta[[1]])))
  expect_equal(a$beta[[1]], ols[[1]] + a$sigma_u * sqrt(2 / pi))
})

test_that("ACOLS refuses the models the paper does not cover", {
  set.seed(4)
  X <- cbind(1, rnorm(60)); y <- rnorm(60)
  for (m in c("NG", "NB")) {
    expect_error(sfa:::.cols_fit(y, X, m, 1L, moment = "absolute"),
      "NHN", info = m)
  }
})

test_that("ACOLS reaches sigma_u = 0 far less often than COLS", {
  ## The paper's headline claim, on its own design: sigma_eps = 1 with
  ## sigma_eps^2 = sigma_v^2 + (pi-2)/pi sigma_u^2, lambda = 0.4, n = 200.
  ## Parmeter and Zhao (2023) Table 1 reports Type I failure of 0.478 for COLS
  ## and 0.146 for ACOLS there; 300 replications is enough to separate those.
  skip_on_cran()
  set.seed(20260911)
  lam <- 0.4
  su <- sqrt(1 / ((pi - 2) / pi + 1 / lam^2)); sv <- su / lam
  n <- 200; R <- 300
  z <- matrix(NA, R, 2)
  for (b in seq_len(R)) {
    x <- rnorm(n); X <- cbind(1, x)
    y <- 1 + x + rnorm(n, 0, sv) - abs(rnorm(n, 0, su))
    z[b, 1] <- sfa:::.cols_fit(y, X, "NHN", 1L, moment = "third")$sigma_u
    z[b, 2] <- sfa:::.cols_fit(y, X, "NHN", 1L, moment = "absolute")$sigma_u
  }
  cols_fail <- mean(z[, 1] <= 0.01)
  acols_fail <- mean(z[, 2] <= 0.01)
  expect_gt(cols_fail, 0.40)            # paper: 0.478
  expect_lt(acols_fail, 0.25)           # paper: 0.146
  expect_lt(acols_fail, cols_fail / 2)
})

test_that("sfm(estimator = \"acols\") wires the branch up", {
  set.seed(5)
  d <- data.frame(x = rnorm(250))
  d$y <- 1 + 0.8 * d$x + rnorm(250, 0, 0.6) - abs(rnorm(250, 0, 1))
  f <- sfm(y ~ x, d, model_name = "NHN", estimator = "acols")
  expect_s3_class(f, "sfareg")
  expect_identical(f$estimator, "acols")
  expect_equal(colnames(f$out), c("par", "st_err", "t-val"))
  expect_equal(rownames(f$out), c("sigv", "sigu", "(Intercept)", "x"))
  expect_true(all(is.finite(f$out[, "par"])))
  ## The analytic Coelli errors are for the THIRD-moment inversion and must
  ## not be reported here.
  expect_true(all(is.na(f$out[c("sigv", "sigu"), "st_err"])))
  ## robust = is a likelihood idea and has to be refused.
  expect_error(sfm(y ~ x, d, model_name = "NHN", estimator = "acols",
    robust = "mdpd"), "moment estimator")
})

## --- constrained MLE --------------------------------------------------------

test_that("the CMLE constraints hold at the fitted values", {
  set.seed(6)
  n <- 400
  x <- rnorm(n)
  X <- cbind("(Intercept)" = 1, x = x)
  y <- 1 + 0.8 * x + rnorm(n, 0, 0.6) - abs(rnorm(n, 0, 1))
  f <- sfa:::.cmle_fit(y, X, "NHN", 1L)
  eps <- as.numeric(y - X %*% f$beta)
  m1 <- mean(abs(eps)); s2 <- mean((eps - mean(eps))^2)
  sg <- sqrt(f$sigma_u^2 + f$sigma_v^2)
  ## Eq. (2.2): E|eps| = sigma sqrt(2/pi)
  expect_equal(m1, sg * sqrt(2 / pi), tolerance = 1e-8)
  ## Eq. (2.3): Var(eps) = sigma_v^2 + (pi-2)/pi sigma_u^2
  expect_equal(s2, f$sigma_v^2 + (pi - 2) / pi * f$sigma_u^2, tolerance = 1e-8)

  y2 <- 1 + 0.8 * x + rnorm(n, 0, 0.6) - rexp(n, 1)
  g <- sfa:::.cmle_fit(y2, X, "NE", 1L)
  e2 <- as.numeric(y2 - X %*% g$beta)
  ## Eq. (2.8) and Proposition 2.1
  expect_equal(mean((e2 - mean(e2))^2), g$sigma_v^2 + g$sigma_u^2, tolerance = 1e-8)
  expect_equal(mean(abs(e2)),
    2 * g$sigma_u * exp(g$sigma_v^2 / (2 * g$sigma_u^2)) *
      pnorm(-g$sigma_v / g$sigma_u) + sqrt(2 / pi) * g$sigma_v,
    tolerance = 1e-7)
})

test_that("CMLE is a CONSTRAINED maximum: never above the unconstrained one", {
  set.seed(7)
  n <- 300
  x <- rnorm(n); X <- cbind(1, x)
  for (mdl in c("NHN", "NE")) {
    u <- if (mdl == "NHN") abs(rnorm(n, 0, 1)) else rexp(n, 1)
    y <- 1 + 0.8 * x + rnorm(n, 0, 0.6) - u
    f <- sfa:::.cmle_fit(y, X, mdl, 1L)
    ## unconstrained ML on the same likelihood closure
    nll <- function(th) {
      v <- sfa:::.cmle_ll(mdl, as.numeric(y - X %*% th[1:2]), exp(th[4]), exp(th[3]))
      if (is.finite(v)) -v else 1e12
    }
    o <- optim(c(f$beta, log(f$sigma_v), log(f$sigma_u)), nll,
      method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-12))
    expect_true(f$loglik <= -o$value + 1e-6, info = mdl)
  }
})

test_that("CMLE survives wrong skew where plain ML collapses", {
  ## Zhao and Parmeter's whole point. Build a sample whose OLS residuals are
  ## positively skewed and check that ML gives sigma_u = 0 while CMLE does not.
  skip_on_cran()
  set.seed(8)
  n <- 120
  found <- FALSE
  for (try in 1:200) {
    x <- rnorm(n); X <- cbind(1, x)
    y <- 1 + x + rnorm(n, 0, 0.95) - abs(rnorm(n, 0, 0.4))
    e <- lm.fit(X, y)$residuals; e <- e - mean(e)
    if (mean(e^3) > 0) { found <- TRUE; break }
  }
  expect_true(found)
  nll <- function(th) {
    v <- sfa:::.cmle_ll("NHN", as.numeric(y - X %*% th[1:2]), exp(th[4]), exp(th[3]))
    if (is.finite(v)) -v else 1e12
  }
  o <- optim(c(lm.fit(X, y)$coefficients, log(0.7), log(0.7)), nll,
    method = "Nelder-Mead", control = list(maxit = 8000, reltol = 1e-12))
  o <- optim(o$par, nll, method = "Nelder-Mead",
    control = list(maxit = 8000, reltol = 1e-12))
  f <- sfa:::.cmle_fit(y, X, "NHN", 1L)
  expect_lt(exp(o$par[4]), 0.01)        # ML on the Waldman boundary
  expect_gt(f$sigma_u, 0.05)            # CMLE is not
})

test_that("sfm(estimator = \"cmle\") wires the branch up", {
  set.seed(9)
  d <- data.frame(x = rnorm(250))
  d$y <- 1 + 0.8 * d$x + rnorm(250, 0, 0.6) - abs(rnorm(250, 0, 1))
  f <- sfm(y ~ x, d, model_name = "NHN", estimator = "cmle")
  expect_s3_class(f, "sfareg")
  expect_identical(f$estimator, "cmle")
  expect_equal(rownames(f$out), c("sigv", "sigu", "(Intercept)", "x"))
  expect_true(all(is.finite(f$out[, "par"])))
  expect_true(all(is.finite(f$out[, "st_err"])))
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_error(sfm(y ~ x, d, model_name = "NG", estimator = "cmle"), "NHN")
})

## --- the third-moment decomposition -----------------------------------------

test_that("the decomposition reproduces the FGM cross moments exactly", {
  ## For the FGM copula c(a,b) = 1 + theta(1-2a)(1-2b) every cross moment
  ## FACTORS, so the 2-D quadrature can be checked against one-dimensional
  ## integrals that share none of its machinery.
  for (p in list(c(0.5, 0.3, 0.6, -0.7), c(1.2, 0.9, 2.5, 0.5))) {
    du <- p[1]; dv <- p[2]; av <- p[3]; th <- p[4]
    A1 <- integrate(function(v) v * (1 - 2 * exp(sfa:::.gl_lp(v, av, dv))) *
      exp(sfa:::.gl_ld(v, av, dv)), -Inf, Inf)$value
    A2 <- integrate(function(v) v^2 * (1 - 2 * exp(sfa:::.gl_lp(v, av, dv))) *
      exp(sfa:::.gl_ld(v, av, dv)), -Inf, Inf)$value
    g <- sfa:::.sd_grid(du, dv, av, "exponential", "glogistic", "fgm", th, 400L)
    cr <- sfa:::.sd_cross(g)
    ## B1 = -du/2 and B2 = -du^2/2 for an exponential u, both exact.
    expect_equal(cr[["v2u"]], th * A2 * (-du / 2), tolerance = 1e-4)
    expect_equal(cr[["vu2"]], th * A1 * (-du^2 / 2), tolerance = 1e-4)
    expect_equal(cr[["vu"]], th * A1 * (-du / 2), tolerance = 1e-4)
  }
})

test_that("the generalized logistic FGM weights have the closed forms we derived", {
  ## A1 = -delta [psi(2a) - psi(a)] and
  ## A2 =  delta^2 [psi'(a) - psi'(2a) - (psi(2a) - psi(a))^2].
  ## These are what the paper's Eq. (10) and (11) should have been built from;
  ## see the note in R/skewness_decomp.R.
  for (av in c(0.4, 1, 2.5, 4)) {
    D <- digamma(2 * av) - digamma(av)
    A1 <- integrate(function(v) v * (1 - 2 * exp(sfa:::.gl_lp(v, av, 1))) *
      exp(sfa:::.gl_ld(v, av, 1)), -Inf, Inf, rel.tol = 1e-12)$value
    A2 <- integrate(function(v) v^2 * (1 - 2 * exp(sfa:::.gl_lp(v, av, 1))) *
      exp(sfa:::.gl_ld(v, av, 1)), -Inf, Inf, rel.tol = 1e-12)$value
    expect_equal(A1, -D, tolerance = 1e-9, info = av)
    expect_equal(A2, trigamma(av) - trigamma(2 * av) - D^2,
      tolerance = 1e-8, info = av)
  }
})

test_that("the three components sum to the composed third moment", {
  ## The corrected closed form, derived in R/skewness_decomp.R:
  ##   m3 = -2 du^3 + dv^3[psi''(a) - psi''(1)]
  ##        + (3/2) th du dv { du D + dv[psi'(a) - psi'(2a) - D^2] }
  m3_exact <- function(du, dv, av, th) {
    D <- digamma(2 * av) - digamma(av)
    -2 * du^3 + dv^3 * (psigamma(av, 2) - psigamma(1, 2)) +
      1.5 * th * du * dv * (du * D + dv * (trigamma(av) - trigamma(2 * av) - D^2))
  }
  for (p in list(c(0.5, 0.3, 0.6, -0.7), c(1.2, 0.9, 2.5, 0.5),
                 c(0.8, 0.6, 1, 0.9), c(1, 0.5, 1.7, 0))) {
    du <- p[1]; dv <- p[2]; av <- p[3]; th <- p[4]
    g <- sfa:::.sd_grid(du, dv, av, "exponential", "glogistic", "fgm", th, 400L)
    cr <- sfa:::.sd_cross(g)
    um <- sfa:::.sd_u_moments(du, "exponential")
    vm <- sfa:::.sd_v_moments(dv, av, "glogistic")
    got <- -um$m3 + vm$m3 + 3 * (cr[["vu2"]] - cr[["v2u"]])
    expect_equal(got, m3_exact(du, dv, av, th), tolerance = 1e-4,
      info = paste(p, collapse = " "))
    ## the variance identity from the same pieces
    expect_equal(vm$var + um$var - 2 * cr[["vu"]],
      du^2 + dv^2 * (trigamma(av) + trigamma(1)) -
        th * du * dv * (digamma(2 * av) - digamma(av)),
      tolerance = 1e-4, info = paste(p, collapse = " "))
  }
})

test_that("the classic model has nowhere to put a positive third moment", {
  ## Bonanno, De Giovanni and Domma's argument in one assertion: under a normal
  ## noise and independence, two of the three components are IDENTICALLY zero,
  ## so the sign of m3 is the sign of -E[(u - Eu)^3] and nothing else.
  set.seed(10)
  d <- data.frame(x = rnorm(300))
  d$y <- 1 + 0.8 * d$x + rnorm(300, 0, 0.6) - abs(rnorm(300, 0, 1))
  f <- sfm(y ~ x, d, model_name = "NHN", estimator = "cols")
  s <- skewness_decomp(f)
  expect_s3_class(s, "sfa_skew_decomp")
  expect_identical(s$components[["noise"]], 0)
  expect_identical(s$components[["dependence"]], 0)
  expect_lt(s$total, 0)
  ## -E[(u - Eu)^3] for a half normal is sigma^3 sqrt(2/pi)(1 - 4/pi), which is
  ## NEGATIVE -- that minus sign is the whole wrong-skewness problem.
  expect_equal(s$components[["inefficiency"]],
    sqrt(2 / pi) * (1 - 4 / pi) * f$out[["sigu", "par"]]^3, tolerance = 1e-10)
  expect_lt(s$components[["inefficiency"]], 0)
  expect_equal(sum(s$components), s$total)
  expect_output(print(s), "Third-moment decomposition")
})

test_that("skewness_decomp refuses what it cannot decompose", {
  set.seed(11)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rnorm(200, 0, 0.6) - abs(rnorm(200, 0, 1))
  f <- sfm(y ~ x, d, model_name = "NG", estimator = "cols")
  expect_error(skewness_decomp(f), "NHN")
  expect_error(skewness_decomp(lm(y ~ x, d)), "sfareg")
})

test_that("skewness_decomp reads BOTH parameterizations sfm() reports", {
  ## NHN under maximum likelihood reports (lambda, sigma); under the moment
  ## and constrained estimators it reports (sigv, sigu). Reading only one of
  ## them is a "subscript out of bounds" on the other, which is how this was
  ## found.
  set.seed(41)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rnorm(200, 0, 0.6) - abs(rnorm(200, 0, 1))
  f_ml <- sfm(y ~ x, d, model_name = "NHN", estimator = "mle")
  expect_true("lambda" %in% rownames(f_ml$out))
  s_ml <- skewness_decomp(f_ml)
  ## the decomposition must agree with the sigmas the reparameterization implies
  lam <- f_ml$out[["lambda", "par"]]
  sg <- f_ml$out[["sigma", "par"]]
  su <- lam * sg / sqrt(1 + lam^2)
  expect_equal(s_ml$components[["inefficiency"]],
    sqrt(2 / pi) * (1 - 4 / pi) * su^3, tolerance = 1e-10)

  f_c <- sfm(y ~ x, d, model_name = "NHN", estimator = "cmle")
  expect_true("sigu" %in% rownames(f_c$out))
  expect_s3_class(skewness_decomp(f_c), "sfa_skew_decomp")
})
