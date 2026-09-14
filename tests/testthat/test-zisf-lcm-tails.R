## Gap A30. zsfm()'s inefficient-regime log-density was log(dnorm()) +
## log(pnorm()), which is -Inf once either factor underflows, and zsfm() and
## lcsfm() computed the JLMS predictor as dnorm(zz) / pnorm(zz), which is
## 0 / 0 = NaN below zz ~ -37. Both are now taken in logs.

test_that(".jlms_u() is finite where the levels ratio was NaN", {
  sig <- 1.05; lambda <- 3.2
  sigv <- sig / sqrt(1 + lambda^2); sigu <- lambda * sigv
  sstar <- sigu * sigv / sig
  mu <- -14 * sigu^2 / sig^2
  z <- mu / sstar
  expect_true(is.nan(mu + sstar * dnorm(z) / pnorm(z)))
  v <- .jlms_u(mu, sstar)
  expect_true(is.finite(v))
  expect_equal(v, mu + sstar * exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE)), tolerance = 1e-12)
})

test_that("zsfm() returns finite efficiencies and posterior probabilities with an extreme residual", {
  skip_on_cran()
  set.seed(8)
  n <- 300; x <- rnorm(n)
  u <- ifelse(runif(n) > 0.4, abs(rnorm(n)), 0)
  y <- 1 + 0.5 * x + rnorm(n, 0, 0.3) - u
  y[1] <- y[1] + 60
  f <- suppressWarnings(zsfm(y ~ x, model_name = "ZISF", data = data.frame(y, x)))
  expect_true(all(is.finite(as.numeric(f$jlms))))
  expect_true(all(is.finite(as.numeric(f$post.prob))))
})
