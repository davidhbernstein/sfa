## Gap A31. zsfm("ZISF_Z") could stop in a spurious local optimum of its regime
## link, 10 to 64 log-likelihood units below the best maximum, in 7 of 50
## replications. It now tries several link intercepts before optimizing.

test_that("ZISF_Z reaches the better regime-link basin on a sample where it used not to", {
  skip_on_cran()
  set.seed(950009)
  n <- 400; x <- rnorm(n); z <- rnorm(n)
  u <- ifelse(runif(n) > 0.4, abs(rnorm(n, 0, 1)), 0)
  y <- 1 + 0.5 * x + rnorm(n, 0, 0.3) - u
  f <- suppressWarnings(zsfm(y ~ x | z, model_name = "ZISF_Z", data = data.frame(y, x, z)))
  ## Truth: 40% efficient (link intercept qlogis(0.4) = -0.41), no dependence on z.
  ## The spurious basin sits near intercept 0.2, slope 0.7, about 43 units lower.
  expect_lt(f$opt$value, 380)
  expect_lt(f$opt$par[5], 0)
  expect_lt(abs(f$opt$par[6]), 0.3)
})
