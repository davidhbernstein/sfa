## Gap A88: efficiency() on a zsfm() fit. zsfm() stored only `jlms`, which its
## help page documents as E[u | e] CONDITIONAL on the inefficient regime, and
## `post.prob`. efficiency(type = "jlms") read that slot as E[u | e] over both
## regimes, so an observation with a 0.97 posterior probability of being fully
## efficient was reported at 0.84; type = "bc" refused outright. The fit now
## also stores the two-regime u_hat and exp_u_hat, which are what the generics
## read. Checked against an independent integral of the mixture posterior AT
## THE FIT'S OWN PARAMETERS, so the assertion is an identity and does not
## depend on where the optimizer lands.

.a88_post <- function(e, p0, sv, su) {
  k <- function(u) dnorm(e + u, 0, sv) * 2 * dnorm(u, 0, su)
  f0 <- dnorm(e, 0, sv)
  a <- integrate(k, 0, Inf, rel.tol = 1e-11)$value
  b <- integrate(function(u) exp(-u) * k(u), 0, Inf, rel.tol = 1e-11)$value
  c <- integrate(function(u) u * k(u), 0, Inf, rel.tol = 1e-11)$value
  den <- p0 * f0 + (1 - p0) * a
  c(bc = (p0 * f0 + (1 - p0) * b) / den, u = (1 - p0) * c / den,
    p0 = p0 * f0 / den, u_cond = c / a)
}

test_that("zsfm(): efficiency() reports the posterior over both regimes (A88)", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 23, sig_u = 1, sig_v = .3, cons = .5,
                   beta1 = .5, beta2 = .5, a = 5, mu = .5)
  f <- suppressWarnings(zsfm(y_zisf ~ x1 + x2, "ZISF", data = d))
  p <- setNames(f$out[, "par"], rownames(f$out))
  p0 <- exp(-abs(p[["gamma"]])); sv <- abs(p[["sigv"]]); su <- abs(p[["sigu"]])
  skip_if(sv < 0.02, "sigma_v on its boundary: the posterior is degenerate")
  e <- as.numeric(residuals(f))

  ## The identities that define the two new slots.
  expect_equal(f$u_hat, as.numeric((1 - f$post.prob) * f$jlms))
  expect_true(all(f$exp_u_hat >= f$post.prob - 1e-12 & f$exp_u_hat <= 1))

  ## Against the integral, at the observations with the most and least
  ## posterior mass on full efficiency, where the two readings differ most.
  idx <- unique(c(order(f$post.prob, decreasing = TRUE)[1:3], order(f$post.prob)[1:2]))
  ref <- t(vapply(e[idx], .a88_post, numeric(4), p0 = p0, sv = sv, su = su))
  expect_equal(f$post.prob[idx], ref[, "p0"], tolerance = 1e-7)
  expect_equal(f$jlms[idx], ref[, "u_cond"], tolerance = 1e-7)  # unchanged meaning
  expect_equal(efficiency(f)[idx], ref[, "bc"], tolerance = 1e-7)
  expect_equal(efficiency(f, type = "jlms")[idx], exp(-ref[, "u"]), tolerance = 1e-7)
  expect_equal(predict(f, type = "efficiency")[idx], ref[, "bc"], tolerance = 1e-7)
})

test_that("zsfm(\"ZISF_Z\"): the two-regime slots follow the per-observation mixing probability", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 23, sig_u = 1, sig_v = .3, cons = .5,
                   beta1 = .5, beta2 = .5, a = 5, mu = .5)
  f <- suppressWarnings(zsfm(y_zisf_z ~ x1 + x2 | z, "ZISF_Z", data = d))
  expect_length(f$exp_u_hat, nrow(d))
  expect_equal(f$u_hat, as.numeric((1 - f$post.prob) * f$jlms))
  expect_equal(efficiency(f, type = "jlms"), exp(-f$u_hat))
})
