## Cornwell-Schmidt-Sickles (1990) and Lee-Schmidt (1993).

.cl_panel <- function(N = 90, Tt = 8, seed = 3) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  n <- nrow(d)
  d$x1 <- runif(n, 1, 5)
  d$x2 <- runif(n, 1, 5)
  d
}

test_that("CSS recovers the frontier under firm-specific quadratic effects", {
  skip_on_cran()
  d <- .cl_panel()
  N <- max(d$name)
  th <- cbind(rnorm(N, 0, .5), rnorm(N, 0, .08), rnorm(N, 0, .008))
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 +
    th[d$name, 1] + th[d$name, 2] * d$t + th[d$name, 3] * d$t^2 +
    rnorm(nrow(d), 0, .2)

  f <- psfm(y ~ x1 + x2, model_name = "CSS", data = d,
    individual = "name", time = "year")

  expect_s3_class(f, "sfareg")
  expect_identical(names(coef(f)), c("x1", "x2"))
  expect_true(all(abs(coef(f) - c(0.6, 0.3)) < 3 * f$std.errors))
  expect_equal(f$sigma_v, 0.2, tolerance = 0.05)
  expect_equal(dim(f$theta), c(N, 3L))
  ## Time-varying, and normalized within each period.
  expect_length(f$u_hat, nrow(d))
  expect_true(all(f$u_hat >= 0))
  expect_true(all(tapply(f$u_hat, d$year, min) < 1e-8))
  ## Not maximum likelihood, so no information criteria.
  expect_true(is.na(suppressWarnings(as.numeric(logLik(f)))))
})

test_that("LS recovers the common temporal pattern and the frontier", {
  skip_on_cran()
  d <- .cl_panel(seed = 9)
  N <- max(d$name)
  delta_true <- c(1, 1.4, 1.9, 2.3, 2.5, 2.4, 2.0, 1.5)
  ai <- -abs(rnorm(N, 0, .6))
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + delta_true[d$t] * ai[d$name] +
    rnorm(nrow(d), 0, .2)

  g <- psfm(y ~ x1 + x2, model_name = "LS", data = d,
    individual = "name", time = "year")

  expect_true(all(abs(coef(g) - c(0.6, 0.3)) < 3 * g$std.errors))
  expect_true(g$ls_converged)
  ## delta_1 is the normalization, so it is 1 exactly.
  expect_equal(unname(g$delta[1]), 1)
  ## The shape is what is identified; correlate rather than compare levels.
  expect_gt(cor(as.numeric(g$delta), delta_true), 0.99)
  expect_gt(abs(cor(g$alpha_i, ai)), 0.95)
  expect_true(all(tapply(g$u_hat, d$year, min) < 1e-8))
})

test_that("CSS refuses what it cannot identify", {
  skip_on_cran()
  d <- .cl_panel(N = 30, Tt = 3)
  d$y <- 0.6 * d$x1 + rnorm(nrow(d), 0, .2)
  ## T = 3 saturates the firm quadratic exactly.
  expect_error(
    psfm(y ~ x1 + x2, "CSS", d, individual = "name", time = "year"),
    "fewer than 4 periods"
  )

  d2 <- .cl_panel(N = 40, Tt = 8)
  d2$trend <- d2$t
  d2$y <- 0.6 * d2$x1 + rnorm(nrow(d2), 0, .2)
  ## A pure time trend is spanned by every firm's own quadratic.
  expect_error(
    psfm(y ~ x1 + trend, "CSS", d2, individual = "name", time = "year"),
    "not separately identified"
  )
})

test_that("both estimators drop the intercept into the firm effect", {
  skip_on_cran()
  d <- .cl_panel(N = 40, Tt = 6)
  d$y <- 2 + 0.6 * d$x1 + 0.3 * d$x2 - abs(rnorm(nrow(d), 0, .4)) +
    rnorm(nrow(d), 0, .2)
  for (m in c("CSS", "LS")) {
    f <- psfm(y ~ x1 + x2, m, d, individual = "name", time = "year")
    expect_false("(Intercept)" %in% names(coef(f)))
    expect_equal(nobs(f), nrow(d))
  }
})
