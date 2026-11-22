## copsfm(vdist = "t"): Student-t noise, gap N5.

test_that("the t marginal's CDF pair is complementary and matches pt()", {
  v <- c(-50, -3, -0.2, 0, 0.2, 3, 50)
  for (nu in c(2.05, 5, 30)) {
    pq <- .cop_v_pq(v, 0.4, nu, "t")
    expect_equal(pq$p, pt(v / 0.4, nu), tolerance = 1e-14)
    expect_equal(pq$q, pt(v / 0.4, nu, lower.tail = FALSE), tolerance = 1e-14)
    expect_equal(.cop_v_ld(v, 0.4, nu, "t"), dt(v / 0.4, nu, log = TRUE) - log(0.4))
    ## the breakpoint quantiles invert the log-odds
    b <- c(-10, -1, 0, 1, 10)
    expect_equal(.cop_v_lo(.cop_v_at_lev(0.4, nu, "t", b), 0.4, nu, "t"), b, tolerance = 1e-8)
  }
})

test_that("the composed density with t noise integrates to 1, with and without a copula", {
  for (cop in c("independent", "gaussian", "clayton")) {
    cp <- switch(cop, independent = NULL, gaussian = 0.4, clayton = 2)
    f <- function(e) exp(.cop_log_dens(e, 1, 1, 0.3, 5, "hnormal", "t", cop, cp))
    expect_equal(integrate(f, -Inf, Inf, rel.tol = 1e-10)$value, 1, tolerance = 1e-7)
  }
})

test_that("t noise approaches normal noise as nu grows", {
  e <- c(-2, -0.7, 0, 0.4)
  ln <- .cop_log_dens(e, 1, 1, 0.3, 1, "hnormal", "normal", "independent", NULL)
  lt <- .cop_log_dens(e, 1, 1, 0.3, 200, "hnormal", "t", "independent", NULL)
  expect_lt(max(abs(lt - ln)), 0.02)
})

test_that("copsfm(vdist = 't') reports nu_v and recovers it", {
  skip_on_cran()
  set.seed(3); n <- 1500
  x <- rnorm(n); u <- abs(rnorm(n)); v <- 0.3 * rt(n, 5)
  d <- data.frame(y = 1 + 0.5 * x + v - u, x = x)
  fit <- copsfm(y ~ x, data = d, copula = "independent", vdist = "t")
  expect_identical(rownames(fit$out), c("(Intercept)", "x", "sigma_u", "delta_v", "nu_v"))
  p <- fit$out[, "par"]
  expect_lt(abs(p[["sigma_u"]] - 1), 0.15)
  expect_lt(abs(p[["delta_v"]] - 0.3), 0.08)
  expect_true(p[["nu_v"]] > 2.05 && p[["nu_v"]] < 20)
  expect_equal(fit$alpha_v, p[["nu_v"]])
  ## t beats normal noise on t data
  fn <- copsfm(y ~ x, data = d, copula = "independent", vdist = "normal")
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(fn)))
  ## skewness_decomp() reads the shape off the fit: symmetric noise, no noise skew
  expect_equal(skewness_decomp(fit)$components[["noise"]], 0)
})
