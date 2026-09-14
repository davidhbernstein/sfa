## Gap A32. ivsfm() computed JLMS as dnorm(zz) / pmax(pnorm(zz), eps). Once
## pnorm(zz) fell below machine epsilon (zz < -8.1, a firm well above the
## frontier) the floor made the ratio vanish, the inefficiency came back
## negative, and the reported efficiency exp(-jlms) exceeded 1.

.iv_tail_data <- function(seed, n = 800) {
  set.seed(seed)
  w1 <- rnorm(n); w2 <- rnorm(n); x1 <- rnorm(n); eta <- rnorm(n)
  v <- 0.5 * (0.6 * eta + sqrt(1 - 0.6^2) * rnorm(n))
  x2 <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
  u <- abs(rnorm(n, 0, 1))
  y <- 0.5 + 0.8 * x1 - 0.6 * x2 + v - u
  ## 20 above the frontier: enough to push C2SLS's zz past -8.1. IVCF and
  ## IVLIML respond by shrinking sigma_u instead and never reach the floor.
  y[1] <- y[1] + u[1] + 20
  data.frame(y = y, x1 = x1, x2 = x2, w1 = w1, w2 = w2)
}

test_that(".jlms_u() matches the log-space JLMS where the floored ratio broke", {
  sst <- 0.4
  zz <- -8.3
  mus <- zz * sst
  expect_lt(mus + sst * dnorm(zz) / pmax(pnorm(zz), .Machine$double.eps), 0)
  expect_equal(.jlms_u(mus, sst), mus + sst * exp(dnorm(zz, log = TRUE) - pnorm(zz, log.p = TRUE)),
    tolerance = 1e-12)
  expect_gt(.jlms_u(mus, sst), 0)
})

test_that("ivsfm() efficiencies stay in [0, 1] with a firm far above the frontier", {
  skip_on_cran()
  d <- .iv_tail_data(11)
  for (mn in c("C2SLS", "IVCF", "IVLIML")) {
    f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = d, model_name = mn))
    expect_true(all(as.numeric(f$jlms) >= 0), info = mn)
    expect_true(all(as.numeric(f$efficiency) <= 1 & as.numeric(f$efficiency) > 0), info = mn)
  }
})
