## A96, A97: psfm("TRE")/("TRE_Z") score E[exp(-u_it) | e_i] under the fitted
## model, integrating over alpha_i ~ N(0, sig_r^2). The helper is pinned against
## a posterior built by integrate() from the model's definition, not against a
## fit; one fit then checks only that psfm() feeds the helper the fit's own
## parameters and residuals.

.tre_by_integrate <- function(e, su, sv, sr) {
  ## firm kernel: phi(a; 0, sr) * prod_t NHN(e_t - a; su_t, sv)
  s <- sqrt(su^2 + sv^2); lam <- su / sv
  lk <- function(a) dnorm(a, 0, sr, log = TRUE) + sum(log(2) - log(s) +
    dnorm((e - a) / s, log = TRUE) + pnorm(-lam * (e - a) / s, log.p = TRUE))
  m <- optimize(lk, c(-20, 20), maximum = TRUE)$maximum; l0 <- lk(m)
  k <- function(a) vapply(a, function(x) exp(lk(x) - l0), 1)
  bc <- function(mu, ss) exp(-mu + ss^2 / 2 + pnorm(mu / ss - ss, log.p = TRUE) -
    pnorm(mu / ss, log.p = TRUE))
  jl <- function(mu, ss) mu + ss * exp(dnorm(mu / ss, log = TRUE) - pnorm(mu / ss, log.p = TRUE))
  den <- integrate(k, -Inf, Inf, rel.tol = 1e-12)$value
  ss <- su * sv / s
  U <- u <- numeric(length(e))
  for (j in seq_along(e)) {
    g <- function(f) function(a) vapply(a, function(x)
      f(-(e[j] - x) * su[j]^2 / s[j]^2, ss[j]) * exp(lk(x) - l0), 1)
    U[j] <- integrate(g(bc), -Inf, Inf, rel.tol = 1e-12)$value / den
    u[j] <- integrate(g(jl), -Inf, Inf, rel.tol = 1e-12)$value / den
  }
  list(U = U, u = u)
}

test_that(".tre_posterior_te equals the integrated TRE posterior (A96, A97)", {
  set.seed(1)
  id <- rep(c("a", "b", "c", "d"), c(5, 1, 8, 3))
  e <- c(rnorm(5, 0.4, 0.5), -1.2, rnorm(8, -0.3, 0.6), c(0.9, 1.1, 0.7))
  ## homoskedastic (TRE) and per-observation sigma_u (TRE_Z)
  for (su in list(0.6, exp(rnorm(length(e), -0.5, 0.4)))) {
    got <- .tre_posterior_te(e, id, su, sig_v = 0.3, sig_r = 0.8)
    su_v <- rep_len(su, length(e))
    for (f in unique(id)) {
      r <- which(id == f)
      ref <- .tre_by_integrate(e[r], su_v[r], 0.3, 0.8)
      expect_equal(got$U[r], ref$U, tolerance = 1e-8)
      expect_equal(got$u[r], ref$u, tolerance = 1e-8)
    }
  }
  ## a firm whose residuals are all positive -- what the old clamp zeroed
  expect_true(all(got$U[id == "d"] < 1))
})

test_that("a long firm and a tight posterior are still integrated accurately", {
  set.seed(2)
  e <- rnorm(60, 1.5, 0.4)       # alpha far from 0, posterior sd ~ 0.05
  got <- .tre_posterior_te(e, rep(1, 60), 0.4, 0.2, sig_r = 3)
  ref <- .tre_by_integrate(e, rep(0.4, 60), 0.2, 3)
  expect_equal(got$U, ref$U, tolerance = 1e-8)
})

test_that("sig_r at zero gives the cross-sectional NHN posterior", {
  e <- c(-0.5, 0.2, 1)
  got <- .tre_posterior_te(e, c(1, 1, 2), 0.6, 0.3, sig_r = 0)
  s2 <- 0.6^2 + 0.3^2
  expect_equal(got$U, .te_battese_coelli(-e * 0.36 / s2, 0.18 / sqrt(s2)))
  ## and the likelihood's sign symmetry in sig_r
  expect_equal(.tre_posterior_te(e, c(1, 1, 2), 0.6, 0.3, -0.7),
               .tre_posterior_te(e, c(1, 1, 2), 0.6, 0.3, 0.7))
})

test_that("psfm('TRE') scores its own residual at its own parameters", {
  skip_on_cran()
  d <- panel_small(t = 5, N = 30)
  f <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = d, model_name = "TRE",
    individual = "name", halton_num = 50))
  p <- f$out[, "par"]
  su <- p[["lambda"]] * p[["sigma"]] / sqrt(1 + p[["lambda"]]^2)
  dd <- f$data
  eps <- dd$y_gtre - p[["(Intercept)"]] - p[["x1"]] * dd$x1 - p[["x2"]] * dd$x2
  te <- .tre_posterior_te(eps, dd$name, su, su / p[["lambda"]], p[["sigr"]])
  ## 1e-7, not 1e-10: the residual rebuilt here differs from psfm()'s in the
  ## last bits, which moves the quadrature's centre, and the 64-node rule is
  ## accurate to ~1e-9 -- that is the arithmetic, not a wiring difference.
  expect_equal(as.numeric(f$U), te$U, tolerance = 1e-7)
  expect_equal(as.numeric(f$u_hat), te$u, tolerance = 1e-7)
  expect_equal(as.numeric(efficiency(f, type = "jlms")), exp(-te$u), tolerance = 1e-7)
})

test_that("GTRE_Z's predictor uses the likelihood's firm-mean sigma_h (A101)", {
  skip_on_cran()
  ## sigma_h,i is the firm MEAN of sqrt(exp(zp'delta_h)) in the likelihood, so
  ## reversing zp's order within each firm changes neither the objective nor,
  ## with the predictor reading the same mean, the scores. Reading the firm's
  ## first row, as the predictor did, makes H move while the fit stays put.
  d <- as.data.frame(data_gen_p(t = 5, N = 40, rand = 3, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  set.seed(5)
  d$zp_gtre <- d$zp_gtre + rnorm(nrow(d), 0, 0.6)     # varies within firm
  d2 <- d
  d2$zp_gtre <- ave(d$zp_gtre, d$name, FUN = rev)
  fit <- function(dd) suppressWarnings(psfm(y_gtre_z ~ x1 + x2 | z_gtre | zp_gtre,
    model_name = "GTRE_Z", data = dd, individual = "name", halton_num = 30, rand.gtre = 7))
  a <- fit(d); b <- fit(d2)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-6)
  expect_equal(as.numeric(a$H), as.numeric(b$H), tolerance = 1e-5)
  expect_equal(as.numeric(a$U), as.numeric(b$U), tolerance = 1e-5)
})
