## Gap A49: NG's log-likelihood at a collapsing sigma_u.
##
## The density carries z^2/4 + log D_{-mu}(z) with z = eps/sigma_v +
## sigma_v/sigma_u, and .log_pcf() returns -z^2/4 plus a term of order log z.
## Added as two numbers, they cancelled. Measured on 200 observations before
## the fix: the summed log-likelihood was off by 0.03 at sigma_u = 1e-7, by +117
## (spuriously better) at 1e-9, and exactly 0 from 1e-10 down -- finite, so the
## non-finite guard let it through, and far above any real optimum.

## log D_nu(z) + z^2/4 from its integral, independently of the series: with
## p = -nu and t = s / z,
##   D_nu(z) e^(z^2/4) = z^(-p) / Gamma(p) int_0^Inf s^(p-1) e^(-s - s^2/(2 z^2)) ds,
## and s = w^(1/p) removes the s^(p-1) singularity at 0 for p < 1.
scaled_ref <- function(p, z) {
  g <- function(w) {
    s <- w^(1 / p)
    exp(-s - s^2 / (2 * z^2))
  }
  -lgamma(p + 1) - p * log(z) +
    log(stats::integrate(g, 0, Inf, rel.tol = 1e-12)$value)
}

test_that(".log_pcf_scaled() matches the integral where the old sum cancelled", {
  for (p in c(0.3, 1, 2.5, 6)) {
    for (z in c(2e3, 1e5, 3e7, 1e10, 1e15)) {
      expect_equal(.log_pcf_scaled(-p, z), scaled_ref(p, z), tolerance = 1e-9,
        info = sprintf("p = %g, z = %g", p, z))
    }
  }
})

test_that(".log_pcf_scaled() is continuous across its switch to the expansion", {
  ## Both forms are valid near the switch; the plain sum loses only about
  ## z^2/4 * 1e-16 there.
  for (p in c(0.5, 2, 8)) {
    zc <- max(1e3, sqrt(1e3) * (p + 2))
    zs <- zc * c(0.999, 1.001)
    expect_equal(.log_pcf_scaled(-p, zs),
      .log_pcf(-p, zs) + zs^2 / 4, tolerance = 1e-9, info = p)
  }
  ## Below the switch it is the old combination exactly.
  z <- c(-3, 0, 0.5, 4, 40)
  expect_identical(.log_pcf_scaled(-1.5, z), .log_pcf(-1.5, z) + z^2 / 4)
})

test_that("NG's likelihood tends to its normal limit as sigma_u collapses", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  f <- sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d, keep_objective = TRUE)
  p <- f$out[, "par"]
  X <- cbind(1, d$x1, d$x2)
  for (su in c(1e-7, 1e-9, 1e-10, 1e-13, 1e-16)) {
    q <- p
    q[2] <- su
    pc <- as.numeric(f$objective(q, per_obs = TRUE))
    ## u ~ Gamma(mu, sigma_u) degenerates to 0, leaving eps ~ N(0, sigma_v^2).
    ## The first omitted term is O(mu sigma_u / sigma_v), below 1e-6 here.
    lim <- dnorm(as.numeric(d$y_pcs_g - X %*% q[4:6]), 0, q[1], log = TRUE)
    expect_false(any(pc == 0), info = su)
    expect_equal(pc, lim, tolerance = 1e-5, info = su)
  }
})

test_that("NG's E[exp(-u) | eps] survives a large z instead of returning NaN", {
  ## exp((z + sigma_v)^2 / 4) / exp(z^2 / 4) was Inf / Inf past z of about 53.
  ## Reference: the posterior u | eps, proportional to
  ## u^(mu - 1) e^(-u / sigma_u) phi((eps + u) / sigma_v), integrated directly.
  sv <- 0.3; su <- 0.006; mu <- 1.7; eps <- 3
  z <- eps / sv + sv / su
  expect_gt(z, 53)
  got <- exp(.log_pcf_scaled(-mu, z + sv) - .log_pcf_scaled(-mu, z))
  lk <- function(u) (mu - 1) * log(u) - u / su + dnorm((eps + u) / sv, log = TRUE)
  m <- max(lk(seq(1e-8, 0.2, length.out = 2e4)))
  num <- stats::integrate(function(u) exp(lk(u) - m - u), 0, Inf, rel.tol = 1e-12)$value
  den <- stats::integrate(function(u) exp(lk(u) - m), 0, Inf, rel.tol = 1e-12)$value
  expect_equal(got, num / den, tolerance = 1e-8)
})
