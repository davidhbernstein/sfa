## npsfm(): nonparametric stochastic frontier models.
## Every test needs np, which is only in Suggests.
##
## Deliberately does NOT call library(np). npsfm() must work with np merely
## installed, not attached -- attaching it here once hid a bug where np's own
## generic gradients() was called unqualified and only resolved because the
## test harness had put np on the search path.

np_data <- function(n = 120, seed = 42, su = 0.6, sv = 0.25) {
  set.seed(seed)
  x1 <- runif(n, 1, 4)
  x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  list(
    d = data.frame(y = m + rnorm(n, 0, sv) - abs(rnorm(n, 0, su)), x1 = x1, x2 = x2),
    m = m, su = su, sv = sv
  )
}

test_that("npsfm() rejects bad input before doing any work", {
  g <- np_data(30)
  ## a pipe formula is an error, not silently ignored
  expect_error(npsfm(y ~ x1 | x2, data = g$d, method = "FLW"), "single-part formula")
  ## SVKZ is half-normal only
  expect_error(npsfm(y ~ x1 + x2, data = g$d, method = "SVKZ", dist = "exp"),
    "normal-half normal"
  )
  expect_error(npsfm(y ~ x1 + x2, data = g$d, method = "nope"), "'arg'")
})

test_that(".flw_neg_cll is finite, positive-lambda only, and minimized near truth", {
  set.seed(3)
  e <- rnorm(400, 0, 0.25) - abs(rnorm(400, 0, 0.6))
  e <- e - mean(e)
  expect_equal(sfa:::.flw_neg_cll(-1, e), sfa:::.SFA_CONSTANTS$MAX_VALUE)
  expect_equal(sfa:::.flw_neg_cll(0, e), sfa:::.SFA_CONSTANTS$MAX_VALUE)
  expect_true(is.finite(sfa:::.flw_neg_cll(2.4, e)))
  best <- optimize(sfa:::.flw_neg_cll, c(1e-4, 50), e = e)$minimum
  ## true lambda is 2.4; the concentrated likelihood should land in the region
  expect_gt(best, 0.8)
  expect_lt(best, 8)
})

test_that("FLW recovers a known normal-half normal frontier", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(250, seed = 11)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn")

  expect_s3_class(f, "npsfareg")
  expect_equal(f$sigma.u, g$su, tolerance = 0.30)
  expect_equal(f$sigma.v, g$sv, tolerance = 0.40)
  expect_equal(f$lambda, g$su / g$sv, tolerance = 0.50)
  ## the frontier tracks the truth
  expect_gt(cor(f$frontier, g$m), 0.85)
  ## and lies above the conditional mean (production frontier)
  expect_true(all(f$frontier >= f$conditional.mean))
  ## efficiency predictions are on (0, 1]
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_true(all(f$u_hat >= 0))
  expect_equal(nobs(f), nrow(g$d))
  expect_length(fitted(f), nrow(g$d))
  expect_equal(residuals(f), g$d$y - f$frontier, ignore_attr = TRUE)
})

test_that("SVKZ returns per-observation sigmas and flags wrong skew", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(250, seed = 11)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "SVKZ")

  expect_length(f$sigma.u, nrow(g$d))
  expect_length(f$sigma.v, nrow(g$d))
  expect_true(all(f$sigma.u >= 0))
  expect_true(all(f$sigma.v >= 0))
  expect_type(f$wrong.skew, "logical")
  ## wrong-skew points are floored at zero and contribute no gradient
  expect_true(all(f$sigma.u[f$wrong.skew] == 0))
  expect_true(all(f$sigma.u.grad[f$wrong.skew, ] == 0))
  expect_equal(dim(f$frontier.grad), c(nrow(g$d), 2L))
})

test_that("cost frontiers recover a cost DGP and sit below the conditional mean", {
  skip_if_not_installed("np")
  skip_on_cran()
  ## A cost frontier adds inefficiency: y = m(x) + v + u.
  set.seed(5); n <- 250
  x1 <- runif(n, 1, 4); x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  dcost <- data.frame(y = m + rnorm(n, 0, 0.25) + abs(rnorm(n, 0, 0.6)),
                      x1 = x1, x2 = x2)
  f <- npsfm(y ~ x1 + x2, data = dcost, method = "FLW", dist = "hn", cost = TRUE)
  expect_equal(f$sigma.u, 0.6, tolerance = 0.40)
  ## the cost frontier is the lower envelope: below the conditional mean
  expect_true(all(f$frontier <= f$conditional.mean))
  expect_gt(cor(f$frontier, m), 0.85)
})

test_that("the wrong orientation warns instead of silently returning sigma_u = 0", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(150, seed = 5)          # production data
  expect_warning(
    f <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn", cost = TRUE),
    "not negatively skewed"
  )
  expect_lt(f$sigma.u, 0.05)
})

test_that("supplying bandwidths skips cross-validation and changes the fit", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(120, seed = 8)
  f1 <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn",
              bw = c(0.5, 0.5))
  f2 <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn",
              bw = c(5, 5))
  expect_false(isTRUE(all.equal(f1$frontier, f2$frontier)))
  expect_equal(as.numeric(f1$bws$bw), c(0.5, 0.5))
})

test_that("the exponential branch inverts its moments", {
  skip_if_not_installed("np")
  skip_on_cran()
  set.seed(21); n <- 250
  x1 <- runif(n, 1, 4); x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  d <- data.frame(y = m + rnorm(n, 0, 0.25) - rexp(n, rate = 1 / 0.6), x1 = x1, x2 = x2)
  f <- npsfm(y ~ x1 + x2, data = d, method = "FLW", dist = "exp")
  expect_equal(f$sigma.u, 0.6, tolerance = 0.45)
  expect_equal(1 / f$theta, f$sigma.u, tolerance = 1e-8)
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
})
