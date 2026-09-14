## Gap A29. NU's log-density took log(Phi(b) - Phi(a)) as a floored difference
## of CDFs, which rounded to 0 for a residual well above the frontier (a
## log-density of -708 where the truth was -44); NHN and NHN_Z floored the
## density in levels, which capped the penalty on extreme trial points.

.ll_i <- function(model, x, d) {
  f <- suppressWarnings(sfm(if (model == "NHN_Z") y ~ x | z else y ~ x,
    model_name = model, data = d, keep_objective = TRUE))
  f$objective(x, per_obs = TRUE)
}

test_that("NU per-observation log-density is exact in both tails", {
  set.seed(3)
  n <- 50; d <- data.frame(x = rnorm(n)); d$y <- 1 + 0.5 * d$x + rnorm(n, 0, 0.2) - runif(n)
  f <- suppressWarnings(sfm(y ~ x, model_name = "NU", data = d, keep_objective = TRUE))
  p <- f$opt$par
  sv <- 0.05; th <- 1
  pp <- p; pp[1] <- sv; pp[2] <- th
  ll <- f$objective(pp, per_obs = TRUE)
  eps <- d$y - pp[3] - pp[4] * d$x
  a <- eps / sv; b <- (eps + th) / sv
  ref <- ifelse(a > 0,
    pnorm(-a, log.p = TRUE) + log1p(-exp(pnorm(-b, log.p = TRUE) - pnorm(-a, log.p = TRUE))),
    pnorm(b, log.p = TRUE) + log1p(-exp(pnorm(a, log.p = TRUE) - pnorm(b, log.p = TRUE)))) - log(th)
  expect_equal(as.numeric(ll), ref, tolerance = 1e-10)
  expect_true(all(ll > -700))
})

test_that("NHN per-observation log-density is not floored at extreme parameters", {
  set.seed(3)
  n <- 50; d <- data.frame(x = rnorm(n)); d$y <- 1 + 0.5 * d$x + rnorm(n, 0, 0.3) - abs(rnorm(n))
  f <- suppressWarnings(sfm(y ~ x, model_name = "NHN", data = d, keep_objective = TRUE))
  pp <- f$opt$par; pp[1] <- 50; pp[2] <- 0.05
  ll <- f$objective(pp, per_obs = TRUE)
  eps <- d$y - pp[3] - pp[4] * d$x
  ref <- log(2) - log(pp[2]) + dnorm(eps / pp[2], log = TRUE) + pnorm(-eps * pp[1] / pp[2], log.p = TRUE)
  expect_equal(as.numeric(ll), ref, tolerance = 1e-10)
  expect_true(any(ll < -800))
})
