## A100: psfm("PL80_MVTN")'s orthant probabilities and scores. The likelihood
## meets only equicorrelated covariances r I + q J; for q >= 0 the orthant is a
## 1-D integral over the common factor, so each piece is pinned against
## integrate() on that integral, not against a fit.

.orth_by_integrate <- function(m, r, q) {
  lw <- function(c) dnorm(c, log = TRUE) +
    vapply(c, function(cc) sum(pnorm(-(m + sqrt(q) * cc) / sqrt(r), log.p = TRUE)), 1)
  md <- optimize(lw, c(-60, 60), maximum = TRUE)$maximum; l0 <- lw(md)
  l0 + log(integrate(function(c) exp(lw(c) - l0), -Inf, Inf, rel.tol = 1e-13)$value)
}

test_that("the equicorrelated orthant is exact for q >= 0, deep in the tail too", {
  set.seed(1)
  cases <- list(list(rnorm(6), 0.5, 0.3), list(rnorm(6, 2), 0.2, 0.5),
                list(rep(4, 6), 0.3, 0.1), list(rnorm(10, -1), 1, 2),
                list(rnorm(3), 0.4, 0))
  for (a in cases) expect_equal(.equi_log_orthant(a[[1]], a[[2]], a[[3]]),
    .orth_by_integrate(a[[1]], a[[2]], a[[3]]), tolerance = 1e-10)
  ## rows are independent problems
  M <- rbind(cases[[1]][[1]], rnorm(6, 1))
  expect_equal(.equi_log_orthant(M, 0.5, 0.3),
    c(.orth_by_integrate(M[1, ], 0.5, 0.3), .orth_by_integrate(M[2, ], 0.5, 0.3)),
    tolerance = 1e-10)
})

test_that("q < 0 (rho < 0) falls back to a converged deterministic rule", {
  set.seed(2); m <- rnorm(6)
  ref <- log(as.numeric(mvtnorm::pmvnorm(upper = rep(0, 6), mean = m,
    sigma = diag(0.5, 6) - 0.05,
    algorithm = mvtnorm::GenzBretz(maxpts = 1e7, abseps = 1e-14, releps = 1e-10))))
  expect_equal(.equi_log_orthant(m, 0.5, -0.05), ref, tolerance = 1e-6)
  ## outside positive-definiteness: NA, never a number
  expect_true(is.na(.equi_log_orthant(m, 0.5, -0.2)))
})

test_that("u_hat and exp_u_hat are the exact truncated-normal posterior moments", {
  r <- 0.3; q <- 0.2; mu <- c(-0.4, 0.3, 0.8, -1, 0.1, 0.5); ss <- sqrt(r)
  sc <- .pl_mvtn_scores(matrix(mu, 1), r, q)
  lw <- function(c) dnorm(c, log = TRUE) +
    vapply(c, function(cc) sum(pnorm(-(mu + sqrt(q) * cc) / ss, log.p = TRUE)), 1)
  den <- integrate(function(c) exp(lw(c)), -Inf, Inf, rel.tol = 1e-13)$value
  ## given c the coordinates are independent N(mu_t + sqrt(q) c, r) on (-Inf, 0]
  mom <- function(f) sapply(seq_along(mu), function(t) integrate(function(c)
    exp(lw(c)) * vapply(c, function(cc) f(mu[t] + sqrt(q) * cc), 1),
    -Inf, Inf, rel.tol = 1e-12)$value / den)
  Ew <- mom(function(m) m - ss * exp(dnorm(m / ss, log = TRUE) - pnorm(-m / ss, log.p = TRUE)))
  Ee <- mom(function(m) exp(m + r / 2 + pnorm((-m - r) / ss, log.p = TRUE) - pnorm(-m / ss, log.p = TRUE)))
  expect_equal(as.numeric(sc$u_hat), -Ew, tolerance = 1e-10)
  expect_equal(as.numeric(sc$exp_u_hat), Ee, tolerance = 1e-10)
  ## and they are different functionals: exp_u_hat is not exp(-u_hat)
  expect_true(all(sc$exp_u_hat > exp(-sc$u_hat)))
})

test_that("psfm('PL80_MVTN') stores the scores at its own parameters", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 4, N = 40, rand = 7, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  f <- suppressWarnings(psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN",
    data = d, individual = "name"))
  p <- f$out[, "par"]
  pr <- .pl_mvtn_parts(p[["sigu"]], p[["rho"]], p[["sigv"]], 4L)
  E <- matrix(d$y_pl_mvtn - p[["(Intercept)"]] - p[["x1"]] * d$x1 - p[["x2"]] * d$x2,
    ncol = 4, byrow = TRUE)
  sc <- .pl_mvtn_scores(E %*% pr$Q / p[["sigv"]]^2, pr$r_Q, pr$q_Q)
  expect_equal(f$u_hat, as.numeric(t(sc$u_hat)), tolerance = 1e-8)
  expect_equal(f$exp_u_hat, as.numeric(t(sc$exp_u_hat)), tolerance = 1e-8)
})
