## Gap A88's two follow-ups. predict(type = "response") used -log(exp_u_hat) as
## u for fits without a stored u_hat -- the Jensen substitution A68 removed from
## efficiency() -- and sfm("NE_Z") stored no efficiency score at all.

.a88_data <- function() {
  set.seed(3)
  n <- 300
  x <- rnorm(n); z <- rnorm(n)
  u <- rexp(n, 1 / exp(-0.7 + 0.5 * z))
  data.frame(y = 1 + 0.5 * x + rnorm(n, 0, 0.3) - u, x = x, z = z)
}

test_that("NE_Z's scores are the posterior moments at the fit's own parameters", {
  d <- .a88_data()
  f <- suppressWarnings(sfm(y ~ x | z, data = d, model_name = "NE_Z"))
  p <- f$out[, "par"]
  sv <- p[["sigma_v"]]
  su <- exp(p[["Zu.(Intercept)"]] + p[["Zu.z"]] * d$z)
  e <- d$y - p[["(Intercept)"]] - p[["x"]] * d$x
  ## Hand-built posterior by quadrature, at a handful of observations.
  for (i in c(1, 50, 123, 300)) {
    k <- function(u) dnorm(e[i] + u, 0, sv) * dexp(u, 1 / su[i])
    den <- integrate(k, 0, Inf, rel.tol = 1e-12)$value
    eu <- integrate(function(u) u * k(u), 0, Inf, rel.tol = 1e-12)$value / den
    te <- integrate(function(u) exp(-u) * k(u), 0, Inf, rel.tol = 1e-12)$value / den
    expect_equal(f$u_hat[i], eu, tolerance = 1e-8)
    expect_equal(f$exp_u_hat[i], te, tolerance = 1e-8)
  }
  expect_equal(efficiency(f), f$exp_u_hat)
  expect_equal(efficiency(f, type = "jlms"), exp(-f$u_hat))
})

test_that("predict(type = 'response') shifts the frontier by E[u | e]", {
  d <- .a88_data()
  for (m in c("NHN", "NE_Z")) {
    fm <- if (m == "NE_Z") y ~ x | z else y ~ x
    f <- suppressWarnings(sfm(fm, data = d, model_name = m))
    u <- -log(efficiency(f, type = "jlms"))
    expect_equal(as.numeric(predict(f, type = "response")), as.numeric(fitted(f) - u))
  }
})

test_that("a fit storing only the Battese-Coelli score is refused, not approximated", {
  d <- .a88_data()
  f <- suppressWarnings(sfm(y ~ x, data = d, model_name = "NHN"))
  f$u_hat <- f$jlms <- f$u_posterior <- NULL
  expect_false(is.null(f$exp_u_hat))
  expect_error(predict(f, type = "response"), "needs E\\[u \\| e\\]")
  expect_error(efficiency(f, type = "jlms"), "jlms")
})
