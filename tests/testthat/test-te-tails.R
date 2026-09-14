## Gap A28. Battese-Coelli efficiency was computed as a ratio of 1 - pnorm()
## terms, which underflowed to 0 / 0 for a residual far above the frontier and
## reported such a firm as fully inefficient; the mean-efficiency closed forms
## overflowed to NaN for large scales.

test_that(".te_battese_coelli() stays accurate where the ratio form underflowed", {
  su <- 1; sv <- 0.1; sig <- sqrt(su^2 + sv^2)
  eps <- c(0.8, 1, 1.5, 5)
  sstar <- su * sv / sig
  te <- .te_battese_coelli(-eps * su^2 / sig^2, rep(sstar, length(eps)))
  inner <- (su / sv) * eps / sig
  ref <- exp(pnorm(sstar + inner, lower.tail = FALSE, log.p = TRUE) -
    pnorm(inner, lower.tail = FALSE, log.p = TRUE) + (su^2 / sig^2) * (eps + 0.5 * sv^2))
  expect_equal(te, pmin(ref, 1), tolerance = 1e-12)
  expect_true(all(te > 0.98))
  old <- (1 - pnorm(sstar + inner)) / pmax(1 - pnorm(inner), .Machine$double.xmin) *
    exp((su^2 / sig^2) * (eps + 0.5 * sv^2))
  expect_true(any(old[-1] == 0))
})

test_that("NHN efficiencies are the log-space Battese-Coelli predictor", {
  set.seed(4)
  n <- 300; x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n, 0, 0.3) - abs(rnorm(n))
  f <- suppressWarnings(sfm(y ~ x, model_name = "NHN", data = data.frame(y, x)))
  p <- f$out[, "par"]
  su <- p[["lambda"]] * p[["sigma"]] / sqrt(1 + p[["lambda"]]^2)
  sv <- su / p[["lambda"]]
  eps <- y - p[["(Intercept)"]] - p[["x"]] * x
  expect_equal(as.numeric(f$exp_u_hat),
    as.numeric(.te_battese_coelli(-eps * su^2 / p[["sigma"]]^2, su * sv / p[["sigma"]])),
    tolerance = 1e-10)
})

test_that("mean-efficiency closed forms are finite at large scales", {
  expect_equal(exp(.log_nr_g(40) - dnorm(40, log = TRUE)),
    integrate(function(u) exp(-u) * (u / 1600) * exp(-u^2 / 3200), 0, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-8)
  expect_equal(exp(log(2) + 40^2 / 2 + pnorm(-40, log.p = TRUE)),
    integrate(function(u) exp(-u) * 2 * dnorm(u, 0, 40), 0, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-8)
})
