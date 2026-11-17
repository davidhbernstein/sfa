## Gap A89: ttsfm("TTHN", vhet = ...) reported its conditional expectations
## (metrics: Eemu.cond, Eeu.cond, Eemw.cond, Eew.cond, and the M-measures built
## from them) from the wrong parameters. metrics.hn() read sigma_u and sigma_w
## at positions fixed for a ONE-slot sigma_v, so with a vhet block sigma_u came
## from a Zv. coefficient. The likelihood indexed correctly, so the estimates
## were right and only the reported posteriors were wrong -- off by up to 1.7.
##
## The identity checked here holds at whatever point the optimizer stops, so a
## short fit is enough: it asserts consistency between the metrics and the
## likelihood's own parameters, not where the optimizer lands.

test_that("TTHN with vhet: metrics are the posterior at the fit's own parameters (A89)", {
  skip_on_cran()
  d <- data_gen_cs(N = 120, rand = 23, sig_u = 1, sig_v = .3, cons = .5,
                   beta1 = .5, beta2 = .5, a = 5, mu = .5)
  f <- suppressWarnings(ttsfm(y_tthn ~ x1 + x2, model_name = "TTHN", data = d,
                              vhet = ~ x1, maxit.bobyqa = 200, maxit.psoptim = 5,
                              maxit.optim = 20, optHessian = FALSE))
  p <- setNames(f$out[, "par"], rownames(f$out))
  sv <- exp(p[["Zv.(Intercept)"]] + p[["Zv.x1"]] * d$x1)
  su <- abs(p[["sigma_u"]]); sw <- abs(p[["sigma_w"]])
  skip_if(min(sv) < 0.02, "sigma_v on its boundary: the posterior is degenerate")
  e <- as.numeric(residuals(f))

  ## e = v + w - u, w ~ |N(0, sw)|, u ~ |N(0, su)|. The density of v + w is a
  ## skew normal, so the posterior of u is a one-dimensional integral.
  post_u <- function(i) {
    s1 <- sqrt(sv[i]^2 + sw^2)
    fvw <- function(m) 2 / s1 * dnorm(m / s1) * pnorm((sw / sv[i]) * m / s1)
    k <- function(u) fvw(e[i] + u) * 2 * dnorm(u, 0, su)
    at <- max(0, -e[i])
    I <- function(g) integrate(g, 0, at, rel.tol = 1e-11)$value +
                     integrate(g, at, Inf, rel.tol = 1e-11)$value
    a <- I(k)
    c(emu = I(function(u) exp(-u) * k(u)) / a)
  }
  idx <- c(which.min(e), which.max(e), which.min(sv), which.max(sv))
  ref <- vapply(idx, post_u, numeric(1))
  expect_equal(as.numeric(f$metrics$Eemu.cond[idx]), unname(ref), tolerance = 1e-6)
  expect_equal(as.numeric(f$metrics$M6), 1 - as.numeric(f$metrics$Eemu.cond))
})
