## A106: selsfm("kts")'s jlms must be the posterior mean of u under the SAME
## phi_i-weighted mixture the likelihood integrates. Pinned against a
## posterior built by hand with integrate(), at arbitrary parameters -- no fit.

test_that(".kts_jlms() is E[u | y, I] under the phi_i-weighted mixture", {
  set.seed(3)
  n <- 6
  x <- rnorm(n); z <- rnorm(n)
  X <- cbind(1, x); Z <- cbind(1, z)
  I <- c(1, 0, 1, 1, 0, 0)
  b0 <- c(0.4, 0.5); b1 <- c(0.9, 0.3)
  sv0 <- 0.4; sv1 <- 0.6; su0 <- 0.5; su1 <- 1.3
  gam <- c(0.2, -0.7); del <- 0.8
  theta <- c(b0, b1, log(c(sv0, sv1, su0, su1)), gam, del)
  mu <- ifelse(I == 1, X %*% b1, X %*% b0)
  y <- as.numeric(mu) - c(0.3, 1.9, 0.05, 2.6, 0.8, -0.4)
  a <- as.numeric(Z %*% gam)
  sv <- ifelse(I == 1, sv1, sv0)

  fhn <- function(u, s) 2 / s * dnorm(u / s)
  A <- function(i, s) integrate(function(u) pnorm(a[i] + del * u) * fhn(u, s), 0, Inf)$value
  want <- vapply(seq_len(n), function(i) {
    phi <- A(i, su0) / (1 + A(i, su0) - A(i, su1))
    k <- function(u, s) {
      pr <- pnorm(a[i] + del * u)
      dnorm(y[i], mu[i] - u, sv[i]) * (if (I[i] == 1) pr else 1 - pr) * fhn(u, s)
    }
    m <- function(s, p) integrate(function(u) u^p * k(u, s), 0, Inf, rel.tol = 1e-10)$value
    (phi * m(su1, 1) + (1 - phi) * m(su0, 1)) / (phi * m(su1, 0) + (1 - phi) * m(su0, 0))
  }, numeric(1))

  got <- sfa:::.kts_jlms(theta, y, X, X, Z, I, sfa:::.kts_nodes(512L), S = 1)
  expect_equal(got, want, tolerance = 1e-6)

  ## The equal-weight version is a different number on this design, so the
  ## test can tell the two apart.
  phis <- vapply(seq_len(n), function(i) A(i, su0) / (1 + A(i, su0) - A(i, su1)), 1)
  expect_gt(max(abs(phis - 0.5)), 0.1)
})
