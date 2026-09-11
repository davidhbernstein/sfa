## Corrected ordinary least squares, sfm(estimator = "cols").

test_that("COLS recovers the DGP for NHN and NE", {
  skip_on_cran()
  d <- data_gen_cs(N = 4000, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)

  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  p <- coef(f)
  expect_equal(unname(p[["sigv"]]), 0.3, tolerance = 0.15)
  expect_equal(unname(p[["sigu"]]), 1.0, tolerance = 0.15)
  expect_equal(unname(p[["(Intercept)"]]), 0.5, tolerance = 0.15)
  expect_equal(unname(p[["x1"]]), 0.5, tolerance = 0.15)

  g <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  q <- coef(g)
  expect_equal(unname(q[["sigv"]]), 0.3, tolerance = 0.2)
  expect_equal(unname(q[["sigu"]]), 1.0, tolerance = 0.2)
})

test_that("the intercept is actually corrected by E[u]", {
  ## The whole content of the "corrected" in COLS. An earlier version located
  ## the intercept column by a name that data_i_vars does not carry (it holds
  ## make.names()-mangled labels), so match() returned NA and the correction
  ## was silently skipped -- leaving a plain OLS intercept that looked
  ## plausible but was E[u] too low.
  skip_on_cran()
  d <- data_gen_cs(N = 3000, rand = 4, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f  <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  b0 <- coef(f)[["(Intercept)"]]
  su <- coef(f)[["sigu"]]
  ols <- stats::lm(y_pcs ~ x1 + x2, data = d)$coefficients[["(Intercept)"]]

  expect_equal(b0 - ols, su*sqrt(2/pi), tolerance = 1e-8)
  expect_gt(b0, ols)                       ## the frontier lies above the mean
  ## and the slopes are left exactly as OLS gave them
  expect_equal(unname(coef(f)[["x1"]]),
               unname(stats::lm(y_pcs ~ x1 + x2, data = d)$coefficients[["x1"]]),
               tolerance = 1e-10)
})

test_that("NE and NHN corrections use their own E[u], not a shared one", {
  ## E[u] = sigma_u*sqrt(2/pi) for the half-normal but sigma_u for the
  ## exponential; using one formula for both would pass the NHN test above and
  ## quietly mis-shift every NE intercept.
  skip_on_cran()
  d <- data_gen_cs(N = 3000, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  g   <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  ols <- stats::lm(y_pcs_e ~ x1 + x2, data = d)$coefficients[["(Intercept)"]]
  expect_equal(coef(g)[["(Intercept)"]] - ols, coef(g)[["sigu"]], tolerance = 1e-8)
})

test_that("wrong-skew residuals are reported, not silently absorbed", {
  ## Olson, Schmidt and Waldman's Type I failure: with positively skewed
  ## residuals the moment equations have no admissible solution.
  set.seed(3)
  n <- 400
  x1 <- stats::runif(n); x2 <- stats::runif(n)
  ## deliberately POSITIVE skew: an "inefficiency" term with the wrong sign
  y  <- 0.5 + 0.5*x1 + 0.5*x2 + stats::rnorm(n, 0, 0.3) + abs(stats::rnorm(n))
  dd <- data.frame(y, x1, x2)

  expect_warning(f <- sfm(y ~ x1 + x2, model_name = "NHN", data = dd,
                          estimator = "cols"), "WRONG")
  expect_equal(unname(coef(f)[["sigu"]]), 0)
  expect_true(f$wrong_skew)
  expect_gte(f$residual_moments[["m3"]], 0)
  ## no efficiency prediction is possible when sigma_u collapses
  expect_true(all(is.na(f$exp_u_hat)))
  ## and the intercept is then left uncorrected, since E[u] = 0
  ols <- stats::lm(y ~ x1 + x2, data = dd)$coefficients[["(Intercept)"]]
  expect_equal(coef(f)[["(Intercept)"]], ols, tolerance = 1e-10)
})

test_that("COLS is deterministic and needs no optimizer", {
  skip_on_cran()
  d <- data_gen_cs(N = 1000, rand = 6, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  a <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  b <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  expect_equal(coef(a), coef(b))
  expect_null(a$opt)                       ## no optimizer output to store
  expect_equal(a$estimator, "cols")
})

test_that("bootstrap standard errors are produced and reproducible", {
  skip_on_cran()
  d <- data_gen_cs(N = 1000, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  a <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols",
           cols_boot = 100, rand.cols = 11)
  b <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols",
           cols_boot = 100, rand.cols = 11)
  expect_equal(a$std.errors, b$std.errors)
  expect_true(all(is.finite(a$std.errors)))
  expect_equal(dim(a$cols_boot_draws), c(100L, length(coef(a))))
  ## Without the bootstrap, NHN now carries Coelli's (1995, Appendix 1)
  ## analytic delta-method errors rather than NA. They were NA until 1.2.0,
  ## which is what this test used to assert.
  c0 <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  expect_true(is.finite(c0$std.errors[["sigv"]]))
  expect_true(is.finite(c0$std.errors[["sigu"]]))
  expect_true(is.finite(c0$std.errors[["(Intercept)"]]))
  expect_true(is.finite(c0$std.errors[["x1"]]))   ## OLS slope SEs are valid
  ## and they agree with the bootstrap, which is the independent check
  expect_equal(c0$std.errors[["sigu"]], a$std.errors[["sigu"]], tolerance = 0.2)
  expect_equal(c0$std.errors[["sigv"]], a$std.errors[["sigv"]], tolerance = 0.2)
  ## The analytic variance is derived for the half-normal only, so the other
  ## COLS models still carry NA and need the bootstrap.
  ce <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  expect_true(is.na(ce$std.errors[["sigu"]]))
})

test_that("the bootstrap restores the caller's RNG stream", {
  skip_on_cran()
  d <- data_gen_cs(N = 400, rand = 8, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  set.seed(99); ref <- stats::rnorm(3)
  set.seed(99)
  invisible(sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
                estimator = "cols", cols_boot = 20, rand.cols = 5))
  expect_equal(stats::rnorm(3), ref)
})

test_that("unsupported models and incompatible options error clearly", {
  skip_on_cran()
  d <- data_gen_cs(N = 300, rand = 9, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NR", data = d,
                   estimator = "cols"), "NHN")
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
                   estimator = "cols", robust = "mlqe"), "moment estimator")
})

test_that("estimator defaults to mle, leaving existing calls untouched", {
  expect_equal(eval(formals(sfm)$estimator)[1], "mle")
})

test_that("estimator = \"mols\" is an exact synonym for \"cols\"", {
  skip_on_cran()
  ## What .cols_fit() implements is Olson, Schmidt and Waldman's MODIFIED OLS,
  ## not Winsten's corrected OLS, so the method is findable under both names.
  ## The two must not drift into being two code paths.
  d <- cs_small(N = 400)
  a <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  b <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "mols")
  expect_equal(a$out, b$out)
  expect_equal(a$exp_u_hat, b$exp_u_hat)
  ## The stored tag stays "cols", so anything reading $estimator downstream --
  ## and every fit saved by an earlier version -- keeps working.
  expect_equal(b$estimator, "cols")
})

test_that("an unknown estimator is rejected rather than silently ignored", {
  expect_error(
    sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = 100),
        estimator = "winsten"),
    "should be one of"
  )
})

## ---------------------------------------------------------------------------
## Carree (2002): a binomial inefficiency (added 2026-09-10, gap L19)
## ---------------------------------------------------------------------------
##
## The one inefficiency distribution in the package that can be skewed either
## way, and the only reason it is here. Every other one-sided law is positively
## skewed, so a production frontier implies a negatively skewed composed error
## and a positive residual skew has nowhere to go but sigma_u = 0. Carree's
## point is that a binomial with p > 1/2 is negatively skewed and reads the
## same sample moment as "most firms are considerably inefficient".

test_that("the binomial inversion returns the moments it was built from", {
  ## Population moments in, (n, p, sigma_v) out -- no sampling error involved,
  ## so any failure is the algebra. Carree Eq. (7):
  ##   m2 = sigma_v^2 + n p (1-p)
  ##   m3 = -n p (1-p) (1 - 2p)
  ##   m4 - 3 m2^2 = n p (1-p) (1 - 6p + 6p^2)
  for (nb in c(6, 20)) {
    for (pb in c(0.25, 0.4, 0.6, 0.75)) {
      for (sv in c(0.8, 2)) {
        vr <- nb * pb * (1 - pb)
        m2 <- sv^2 + vr
        m3 <- -vr * (1 - 2 * pb)
        k4 <- vr * (1 - 6 * pb + 6 * pb^2)
        ## Invert exactly as .cols_fit() does, from the moments alone.
        xr <- k4 / m3
        rt <- sqrt(xr^2 + 3) / 6
        pp <- if (xr < -1) 0.5 + xr / 6 + rt else if (xr > 1) 0.5 + xr / 6 - rt else
          if (m3 > 0) 0.5 + xr / 6 + rt else 0.5 + xr / 6 - rt
        nn <- -m3 / (pp * (1 - pp) * (1 - 2 * pp))
        expect_equal(pp, pb, tolerance = 1e-8, info = paste(nb, pb, sv))
        expect_equal(nn, nb, tolerance = 1e-6, info = paste(nb, pb, sv))
        expect_equal(sqrt(m2 - nn * pp * (1 - pp)), sv, tolerance = 1e-6)
      }
    }
  }
})

test_that("COLS recovers a binomial inefficiency from data", {
  skip_on_cran()
  set.seed(3)
  n <- 3e5; nb <- 20; pb <- 0.75; sv <- 1.5
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n, 0, sv) - rbinom(n, nb, pb)
  d <- data.frame(y = y, x1 = x1, x2 = x2)
  f <- sfm(y ~ x1 + x2, data = d, model_name = "NB", estimator = "cols")
  p <- coef(f)

  ## p = 0.75 > 1/2, so this sample is skewed the "wrong" way for every other
  ## distribution in the package -- and it is the ADMISSIBLE case here.
  expect_gt(f$residual_moments[["m3"]], 0)
  expect_false(f$wrong_skew)
  expect_equal(unname(p[["p_bin"]]), pb, tolerance = 0.06)
  expect_equal(unname(p[["n_bin"]]), nb, tolerance = 3)
  expect_equal(unname(p[["sigv"]]), sv, tolerance = 0.3)
  expect_equal(unname(p[["sigu"]]), sqrt(nb * pb * (1 - pb)), tolerance = 0.3)
  ## The slopes are OLS, untouched, and the intercept is corrected by E[u] = np.
  expect_equal(unname(p[["x1"]]), 0.8, tolerance = 0.02)
  expect_equal(unname(p[["(Intercept)"]]) -
      stats::lm(y ~ x1 + x2, data = d)$coefficients[["(Intercept)"]],
    unname(p[["n_bin"]] * p[["p_bin"]]), tolerance = 1e-8)
  ## Four parameters beyond the frontier, so se_v must carry four slots. A
  ## single NA for "the extras" silently shortened it while NG, with one extra,
  ## was the only model that had any.
  expect_identical(rownames(f$out)[1:4], c("sigv", "sigu", "n_bin", "p_bin"))
  expect_equal(nrow(f$out), 7L)
})

test_that("the infeasible region is reported as such, not as an estimate", {
  ## Carree: when the excess fourth moment exceeds the third in absolute value
  ## the implied n is negative and there is no binomial solution. HEAVY-TAILED
  ## noise puts a design there reliably, because k4 grows without m3 following.
  skip_on_cran()
  set.seed(1)
  n <- 1500
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 0.5 + 0.8 * d$x1 - 0.4 * d$x2 + rt(n, df = 3) - abs(rnorm(n, 0, 0.3))
  expect_warning(
    f <- sfm(y ~ x1 + x2, data = d, model_name = "NB", estimator = "cols"),
    "no binomial solution"
  )
  expect_true(f$wrong_skew)
  expect_true(is.na(coef(f)[["n_bin"]]))
  expect_equal(unname(coef(f)[["sigu"]]), 0)
  expect_true(all(is.na(f$exp_u_hat)))
  ## The moments really are in Carree's region, not merely somewhere odd.
  m <- f$residual_moments
  expect_gt(m[["m4"]] - 3 * m[["m2"]]^2, abs(m[["m3"]]))
})

test_that("pure noise does NOT make the binomial refuse, which is the catch", {
  ## Worth pinning because the intuition points the other way. With no
  ## inefficiency at all, m3 and k4 are both sampling noise about zero, their
  ## RATIO is an arbitrary number, and the inversion will usually succeed and
  ## return a confident-looking (n, p). The estimator has no way to say "there
  ## is nothing here" -- unlike every other distribution, where a positive m3
  ## is itself the signal. Read a binomial fit alongside a wrong-skewness test
  ## (see skewness_test), not on its own.
  skip_on_cran()
  set.seed(11)
  n <- 2000
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 0.5 + 0.8 * d$x1 - 0.4 * d$x2 + rnorm(n)
  f <- sfm(y ~ x1 + x2, data = d, model_name = "NB", estimator = "cols")
  expect_false(f$wrong_skew)
  expect_false(is.na(coef(f)[["n_bin"]]))
  expect_gt(unname(coef(f)[["sigu"]]), 0)
})

test_that("NB refuses maximum likelihood rather than failing inside it", {
  d <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  d$y <- d$x1 + rnorm(50)
  expect_error(sfm(y ~ x1 + x2, data = d, model_name = "NB"),
    "estimated by corrected OLS")
})

test_that("the COLS efficiency predictor is the one that model implies", {
  ## REGRESSION. Every estimator = "cols" fit used the normal/half-normal
  ## posterior for E[exp(-u)|eps] regardless of model_name, so an NE fit was
  ## given the half-normal answer. The posteriors differ: NHN is truncated
  ## normal at mu* = -eps s2u/(s2u+s2v) with sd sigma_u sigma_v/sigma, NE at
  ## mu* = -eps - sigma_v^2/sigma_u with sd sigma_v.
  skip_on_cran()
  d <- data_gen_cs(N = 3000, rand = 9, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  su <- coef(f)[["sigu"]]; sv <- coef(f)[["sigv"]]
  eps <- as.numeric(d$y_pcs_e - cbind(1, d$x1, d$x2) %*%
    coef(f)[c("(Intercept)", "x1", "x2")])
  want <- sfa:::.te_battese_coelli(-eps - sv^2 / su, rep_len(sv, length(eps)))
  expect_equal(f$exp_u_hat, want, tolerance = 1e-10)
  ## and it is NOT the half-normal answer, which is what it used to return.
  s2u <- su^2; s2v <- sv^2
  nhn <- sfa:::.te_battese_coelli(-eps * s2u / (s2u + s2v),
    su * sv / sqrt(s2u + s2v))
  expect_gt(max(abs(f$exp_u_hat - nhn)), 1e-3)

  ## NG goes through the parabolic-cylinder form, and must reach it: `lnDv` is
  ## a local defined inside the MAXIMUM-LIKELIHOOD path, which the COLS branch
  ## returns long before, so writing the ML expression verbatim here errors
  ## with "could not find function".
  g <- sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d, estimator = "cols")
  pg <- coef(g)
  epsg <- as.numeric(d$y_pcs_g - cbind(1, d$x1, d$x2) %*%
    pg[c("(Intercept)", "x1", "x2")])
  zg <- epsg / pg[["sigv"]] + pg[["sigv"]] / pg[["sigu"]]
  wantg <- exp(((zg + pg[["sigv"]]) / 2)^2 - (zg / 2)^2 +
    sfa:::.log_pcf(-pg[["mu"]], zg + pg[["sigv"]]) - sfa:::.log_pcf(-pg[["mu"]], zg))
  expect_equal(g$exp_u_hat, pmin(pmax(wantg, 0), 1), tolerance = 1e-12)
  expect_true(all(is.finite(g$exp_u_hat)))
})
