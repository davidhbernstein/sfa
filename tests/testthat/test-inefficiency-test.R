## Coelli (1995) Section 3 and Appendix 1.

.ct_data <- function(N, seed, sig_u = 1, sig_v = 0.5) {
  as.data.frame(data_gen_cs(N = N, rand = seed, sig_u = sig_u, sig_v = sig_v,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1
  ))
}

test_that("the chi-bar-square mixture is the documented one", {
  ## One boundary restriction: half the chi2(1) tail, because chi2(0) is a
  ## point mass at zero and contributes nothing above it.
  expect_equal(.chi_bar_p(3.0, 1L), 0.5 * pchisq(3.0, 1, lower.tail = FALSE))
  ## Coelli's own statement of the rule: the size-alpha critical value equals
  ## the chi2(1) critical value for size 2*alpha, i.e. 2.71 at 5%.
  cv <- qchisq(2 * 0.05, df = 1, lower.tail = FALSE)
  expect_equal(cv, 2.705543, tolerance = 1e-5)
  expect_equal(.chi_bar_p(cv, 1L), 0.05, tolerance = 1e-8)
  ## Two restrictions: 0.25/0.5/0.25 over chi2(0,1,2).
  expect_equal(
    .chi_bar_p(3.0, 2L),
    0.5 * pchisq(3, 1, lower.tail = FALSE) + 0.25 * pchisq(3, 2, lower.tail = FALSE)
  )
  expect_equal(.chi_bar_p(0, 1L), 1)
  expect_equal(.chi_bar_p(-1, 1L), 1)
})

test_that("the one-sided LR p-value is exactly half the naive one", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = .ct_data(300, 7))
  tt <- inefficiency_test(f, test = c("lr_1sided", "lr"))
  expect_equal(tt$statistic[1], tt$statistic[2])
  expect_equal(tt$p.value[1], 0.5 * tt$p.value[2], tolerance = 1e-10)
  ## and the naive test can therefore never reject when the one-sided does not
  expect_true(tt$p.value[1] <= tt$p.value[2])
})

test_that("inefficiency_test() returns a well-formed table and rejects on strong data", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = .ct_data(600, 11, sig_u = 2, sig_v = 0.5))
  tt <- inefficiency_test(f)
  expect_s3_class(tt, "data.frame")
  expect_equal(nrow(tt), 4L)
  expect_named(tt, c("test", "statistic", "null", "p.value", "reject"))
  expect_true(all(tt$p.value >= 0 & tt$p.value <= 1))
  ## sigma_u = 4 sigma_v is not a subtle alternative
  expect_true(tt$reject[tt$test == "LR (one-sided)"])
  expect_true(is.finite(attr(tt, "logLik_H1")))
  expect_true(attr(tt, "logLik_H1") >= attr(tt, "logLik_H0") - 1e-6)
})

test_that("inefficiency_test() refuses models Coelli's mixture is not stated for", {
  skip_on_cran()
  f <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = .ct_data(200, 3))
  expect_error(inefficiency_test(f), "half-normal")
})

test_that("cols_sfm() reproduces Coelli's A13 inversion", {
  ## A13 written out directly from the paper, against the package's own
  ## moment inversion -- they are the same estimator.
  set.seed(4)
  e <- rnorm(50000, 0, 0.6) - abs(rnorm(50000, 0, 1.2))
  e <- e - mean(e)
  m2 <- mean(e^2); m3 <- mean(e^3)
  su2_paper <- (m3 / (sqrt(2 / pi) * (1 - 4 / pi)))^(2 / 3)
  sig2_paper <- m2 + (2 / pi) * su2_paper
  ts <- .gtre_two_step(e, e, 0)
  expect_equal(ts$sigmaSq_uv, sig2_paper, tolerance = 1e-10)
  expect_equal(ts$gamma_uv, su2_paper / sig2_paper, tolerance = 1e-10)
})

test_that("cols_sfm() is consistent and its standard errors are not the OLS ones", {
  skip_on_cran()
  R <- 60
  g <- su <- b0 <- numeric(R)
  for (r in seq_len(R)) {
    cc <- cols_sfm(y_pcs ~ x1 + x2, data = .ct_data(1500, 4000 + r))
    g[r] <- cc$gamma; su[r] <- cc$sigma_u; b0[r] <- cc$coefficients[["(Intercept)"]]
  }
  expect_equal(mean(g), 0.8, tolerance = 0.05)    # sigma_u=1, sigma_v=0.5
  expect_equal(mean(su), 1.0, tolerance = 0.08)
  expect_equal(mean(b0), 0.5, tolerance = 0.08)   # intercept shifted by E[u]

  cc <- cols_sfm(y_pcs ~ x1 + x2, data = .ct_data(400, 7))
  expect_true(is.finite(cc$gamma_se) && cc$gamma_se > 0)
  expect_true(is.finite(cc$sigmaSq_se) && cc$sigmaSq_se > 0)
  ## COLS leaves the slopes alone, so only the intercept moves.
  ols <- lm(y_pcs ~ x1 + x2, data = .ct_data(400, 7))
  expect_equal(unname(cc$coefficients[-1]), unname(coef(ols)[-1]), tolerance = 1e-10)
  expect_gt(cc$coefficients[["(Intercept)"]], coef(ols)[["(Intercept)"]])
})

test_that("a cost frontier shifts the intercept the other way", {
  skip_on_cran()
  ## Build COST data: y = xb + v + u, so the composed error is POSITIVELY
  ## skewed and the frontier lies below the mean.
  set.seed(21)
  n <- 3000
  dc <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dc$y <- 0.5 + 0.5 * dc$x1 + 0.5 * dc$x2 + rnorm(n, 0, 0.5) + abs(rnorm(n, 0, 1))
  cost <- cols_sfm(y ~ x1 + x2, data = dc, inefdec = FALSE)
  ols_b0 <- coef(lm(y ~ x1 + x2, data = dc))[["(Intercept)"]]
  expect_lt(cost$coefficients[["(Intercept)"]], ols_b0)
  expect_false(cost$wrong_skew)
  expect_equal(cost$sigma_u, 1, tolerance = 0.15)
  expect_equal(cost$coefficients[["(Intercept)"]], 0.5, tolerance = 0.15)

  ## Reading the SAME data as a production frontier is the Type I failure, and
  ## it must say so rather than returning a silent sigma_u = 0.
  expect_warning(bad <- cols_sfm(y ~ x1 + x2, data = dc, inefdec = TRUE),
                 "wrong way|WRONG WAY")
  expect_true(bad$wrong_skew)
  expect_equal(bad$sigma_u, 0)
})
