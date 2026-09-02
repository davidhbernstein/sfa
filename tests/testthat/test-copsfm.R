## Copula stochastic frontier: copsfm(), gap K5.

cop_gen <- function(seed, n, rho, su = 1, sv = 0.4, b = c(0.5, 0.8, -0.4)) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ## Gaussian copula between the marginals: draw correlated normals, push each
  ## through its own marginal quantile function. w1 = Phi(z1) is F_V(v) and
  ## w2 = Phi(z2) is F_U(u), so (w1, w2) carry exactly the dependence rho.
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  v <- sv * z1
  u <- su * qnorm((1 + pnorm(z2)) / 2)
  data.frame(y = b[1] + b[2] * x1 + b[3] * x2 + v - u, x1 = x1, x2 = x2)
}

test_that("each copula density integrates to 1 and is 1 at independence", {
  ## The check that catches a mistranscribed density. Both are asserted for
  ## every family, which is why families whose density could not be verified
  ## are deliberately absent from .cop_logc().
  gl <- .gauss_legendre_01(80L)
  G <- expand.grid(a = gl$nodes, b = gl$nodes)
  W <- as.vector(outer(gl$weights, gl$weights))
  for (fam in c("gaussian", "fgm")) {
    for (p in c(-0.6, -0.2, 0, 0.3, 0.7)) {
      I <- sum(W * exp(.cop_logc(G$a, G$b, p, fam)))
      expect_equal(I, 1, tolerance = 1e-4, info = paste(fam, p))
    }
    ## At the independence parameter the density is exactly 1, so log c is 0.
    expect_equal(max(abs(.cop_logc(G$a, G$b, 0, fam))), 0, info = fam)
  }
})

test_that("the composed density reduces to normal/half-normal at independence", {
  ## With c = 1 the quadrature must reproduce the closed form the rest of the
  ## package uses. This checks the composition and the change of variable
  ## u = t/(1-t) together.
  gl <- .gauss_legendre_01(128L)
  tt <- gl$nodes; wt <- gl$weights
  uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  su <- 1; sv <- 0.4; e <- seq(-4, 2, length.out = 9)

  U <- matrix(uu, length(e), length(uu), byrow = TRUE)
  V <- e + U
  lg <- dnorm(V, 0, sv, log = TRUE) + log(2) + dnorm(U, 0, su, log = TRUE)
  lg <- sweep(lg, 2, log(wt) + log(jac), "+")
  got <- .log_row_sum_exp(lg)

  sig <- sqrt(su^2 + sv^2); lam <- su / sv
  want <- log(2) - log(sig) + dnorm(e / sig, log = TRUE) +
    pnorm(-e * lam / sig, log.p = TRUE)
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("copsfm returns a well-formed sfareg", {
  skip_on_cran()
  d <- cop_gen(4, 600, 0.5)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "COP")
  expect_identical(f$copula, "gaussian")
  expect_equal(ncol(f$out), 3L)
  expect_identical(rownames(f$out),
    c("(Intercept)", "x1", "x2", "sigma_u", "sigma_v", "rho"))
  expect_gt(f$out["sigma_u", "par"], 0)
  expect_lt(abs(f$out["rho", "par"]), 1)
  expect_equal(length(f$jlms), nrow(d))
  expect_equal(f$efficiency, exp(-f$jlms))
  expect_equal(nobs(f), nrow(d))
  ## fgm names its parameter differently, and reports it.
  g <- copsfm(y ~ x1 + x2, data = d, copula = "fgm")
  expect_true("theta" %in% rownames(g$out))
})

test_that("the frontier slopes are recovered even where rho is not", {
  skip_on_cran()
  ## The dependence parameter needs a large sample (see ?copsfm); the SLOPES do
  ## not, and separating the two is the point of this test.
  d <- cop_gen(4, 2000, 0.6)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.08)
  expect_equal(unname(f$out["x2", "par"]), -0.4, tolerance = 0.08)
})

test_that("copsfm rejects malformed calls", {
  d <- cop_gen(2, 300, 0.3)
  expect_error(copsfm(~ x1 + x2, data = d), "two-sided formula")
  expect_error(copsfm(y ~ x1 | x2, data = d), "must not contain a `\\|` segment")
  expect_error(copsfm(y ~ x1 + x2, data = as.matrix(d)), "must be a data.frame")
  expect_error(copsfm(y ~ x1 + x2, data = d, inefdec = NA), "must be TRUE or FALSE")
  expect_error(copsfm(y ~ x1 + x2, data = d, n_nodes = 4), ">= 16")
  expect_error(copsfm(y ~ x1 + x2, data = d, copula = "clayton"))
})
