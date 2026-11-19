## A105: selsfm("greene") scores the exact 1-D posterior of u at the fitted
## parameters, by quadrature, not by weighting the likelihood's prior draws.
## Pinned against integrate() on the log-scale posterior kernel.

.greene_by_integrate <- function(r0, a, su, sv, rho, S) {
  lk <- function(u) dnorm(u / su, log = TRUE) + dnorm((r0 + S * u) / sv, log = TRUE) +
    pnorm((rho * (r0 + S * u) / sv + a) / sqrt(1 - rho^2), log.p = TRUE)
  sc <- su * sv / sqrt(su^2 + sv^2)
  md <- optimize(lk, c(0, 50), maximum = TRUE)$maximum; l0 <- lk(md)
  lo <- max(0, md - 40 * sc); hi <- md + 40 * sc
  I <- function(f) integrate(function(u) f(u) * exp(lk(u) - l0), lo, hi,
    rel.tol = 1e-12, subdivisions = 2000)$value
  k0 <- I(function(u) 1)
  c(u = I(identity) / k0, bc = I(function(u) exp(-u)) / k0)
}

test_that("the greene posterior moments are exact, including near |rho| = 1", {
  g <- expand.grid(r0 = c(-3, -0.5, 0, 0.7, 2), a = c(-1.5, 0.3, 2),
                   rho = c(-0.99, -0.5, 0, 0.6, 0.99), S = c(1, -1))
  for (rho in unique(g$rho)) for (S in c(1, -1)) {
    G <- g[g$rho == rho & g$S == S, ]
    got <- .selsfm_greene_posterior(G$r0, G$a, 0.8, 0.4, rho, S)   # vectorised, per-row a_i
    ref <- t(mapply(.greene_by_integrate, G$r0, G$a, MoreArgs = list(su = 0.8, sv = 0.4, rho = rho, S = S)))
    expect_equal(unname(got$u), unname(ref[, "u"]), tolerance = 1e-8)
    expect_equal(unname(got$bc), unname(ref[, "bc"]), tolerance = 1e-8)
  }
  ## a posterior of width ~1e-3 pressed against u = 0 (the case a window sized
  ## by s* alone missed by a factor of two)
  expect_equal(unname(.selsfm_greene_posterior(-3, -1.5, 0.8, 0.4, 0.99, -1)$u),
    unname(.greene_by_integrate(-3, -1.5, 0.8, 0.4, 0.99, -1)[["u"]]), tolerance = 1e-8)
})

test_that("selsfm('greene') stores those moments at its own parameters", {
  skip_on_cran()
  set.seed(11)
  n <- 400
  z1 <- rnorm(n); x1 <- rnorm(n)
  v <- rnorm(n, 0, 0.4); w <- 0.6 * v / 0.4 + sqrt(1 - 0.36) * rnorm(n)
  I <- as.integer(0.3 + z1 + w > 0)
  y <- 1 + 0.5 * x1 + v - abs(rnorm(n, 0, 0.8))
  d <- data.frame(y, x1, z1, I)
  f <- suppressWarnings(selsfm(selection = I ~ z1, frontier = y ~ x1, data = d,
    model_name = "greene", Nsim = 200))
  p <- f$out[, "par"]
  ds <- d[d$I == 1, ]
  r0 <- ds$y - p[["(Intercept)"]] - p[["x1"]] * ds$x1
  a_i <- as.numeric(cbind(1, ds$z1) %*% f$alpha)
  post <- .selsfm_greene_posterior(r0, a_i, abs(p[["sigma_u"]]), abs(p[["sigma_v"]]), p[["rho"]], f$S)
  expect_equal(f$jlms, post$u, tolerance = 1e-10)
  expect_equal(f$exp_u_hat, post$bc, tolerance = 1e-10)
  expect_equal(f$efficiency, exp(-f$jlms))
})
