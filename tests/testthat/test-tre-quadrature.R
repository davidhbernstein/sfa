## psfm(integration = "quadrature"): the exact TRE / TRE_Z likelihood by
## adaptive Gauss-Hermite over the firm effect.

## Independent reference: integrate() per firm, centred on the mode.
tre_ref <- function(e, gid, lam, sig, sig_r) {
  n <- length(e)
  lam <- rep_len(lam, n); sig <- rep_len(sig, n)
  vapply(seq_len(max(gid)), function(i) {
    k <- gid == i
    lf <- function(d) vapply(d, function(dd) {
      z1 <- (e[k] - sig_r * dd) / sig[k]
      sum(log(2 / sig[k]) + dnorm(z1, log = TRUE) + pnorm(-lam[k] * z1, log.p = TRUE))
    }, 1) + dnorm(d, log = TRUE)
    m <- optimize(lf, c(-40, 40), maximum = TRUE, tol = 1e-12)
    m$objective + log(integrate(function(d) exp(lf(d) - m$objective),
                                m$maximum - 30, m$maximum + 30,
                                rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value)
  }, 1)
}

test_that(".tre_quad_loglik() is the TRE integral, scalar and per-row scales", {
  set.seed(3)
  N <- 30; Tt <- 6
  gid <- rep(seq_len(N), each = Tt)
  e <- rnorm(N * Tt, sd = 0.4) - abs(rnorm(N * Tt)) + rep(rnorm(N, sd = 0.5), each = Tt)
  gh <- sfa:::.gauss_hermite_nodes(24L)
  expect_equal(sfa:::.tre_quad_loglik(e, gid, N, 2.5, 1.1, 0.5, gh),
               tre_ref(e, gid, 2.5, 1.1, 0.5), tolerance = 1e-10)
  su <- exp(rnorm(N * Tt, sd = 0.3))
  expect_equal(sfa:::.tre_quad_loglik(e, gid, N, su / 0.4, sqrt(0.16 + su^2), 0.7, gh),
               tre_ref(e, gid, su / 0.4, sqrt(0.16 + su^2), 0.7), tolerance = 1e-10)
  ## sigma_r = 0: the firm effect vanishes and the product is exact.
  nhn <- rowsum(log(2 / 1.1) + dnorm(e / 1.1, log = TRUE) + pnorm(-2.5 * e / 1.1, log.p = TRUE),
                gid)[, 1]
  expect_equal(sfa:::.tre_quad_loglik(e, gid, N, 2.5, 1.1, 0, gh), unname(nhn),
               tolerance = 1e-12)
})

test_that("unbalanced firms are integrated over their own periods", {
  set.seed(4)
  gid <- rep(1:5, times = c(1, 3, 7, 2, 12))
  e <- rnorm(length(gid)) - abs(rnorm(length(gid)))
  gh <- sfa:::.gauss_hermite_nodes(24L)
  expect_equal(sfa:::.tre_quad_loglik(e, gid, 5, 1.8, 0.9, 0.6, gh),
               tre_ref(e, gid, 1.8, 0.9, 0.6), tolerance = 1e-10)
})

pd_tre <- function(N = 40, t = 6, rand = 7) {
  d <- data_gen_p(t = t, N = N, rand = rand, sig_u = 1, sig_v = 0.3, sig_r = 0.4,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  plm::pdata.frame(d, index = c("name", "year"))
}

test_that("psfm('TRE', quadrature) evaluates the exact likelihood, both orientations", {
  skip_on_cran()
  pd <- pd_tre()
  for (prod in c(TRUE, FALSE)) {
    p <- pd
    if (!prod) p$y_tre <- -p$y_tre
    f <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", p, individual = "name",
                               integration = "quadrature", inefdec = prod,
                               keep_objective = TRUE))
    expect_identical(f$integration, "quadrature")
    x <- f$opt$par
    e <- (if (prod) 1 else -1) * (as.numeric(p$y_tre) -
      as.numeric(cbind(1, p$x1, p$x2) %*% x[4:6]))
    gid <- as.integer(factor(as.character(p$name), levels = unique(as.character(p$name))))
    expect_equal(f$objective(x, per_obs = TRUE), tre_ref(e, gid, x[1], x[2], x[3]),
                 tolerance = 1e-10, info = paste("prod =", prod))
    expect_equal(as.numeric(logLik(f)), sum(tre_ref(e, gid, x[1], x[2], x[3])),
                 tolerance = 1e-10)
  }
})

test_that("a cost fit of -y is the mirror of the production fit of y", {
  skip_on_cran()
  ## Exact integration over a symmetric alpha has no draw-set asymmetry, so
  ## the two objectives are identical functions of the parameters.
  pd <- pd_tre(N = 30, t = 5, rand = 11)
  a <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pd, individual = "name",
                             integration = "quadrature", keep_objective = TRUE,
                             maxit.bobyqa = 30, maxit.optim = 1, optHessian = FALSE))
  pc <- pd; pc$y_tre <- -pc$y_tre
  b <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pc, individual = "name",
                             integration = "quadrature", inefdec = FALSE,
                             keep_objective = TRUE, maxit.bobyqa = 30,
                             maxit.optim = 1, optHessian = FALSE))
  x <- a$opt$par
  xm <- x; xm[4:6] <- -xm[4:6]
  expect_equal(b$objective(xm), a$objective(x), tolerance = 1e-12)
})

test_that("TRE_Z quadrature reduces to TRE when sigma_u is constant", {
  skip_on_cran()
  pd <- pd_tre(N = 30, t = 5, rand = 5)
  pd$one <- 1
  z <- suppressWarnings(psfm(y_tre ~ x1 + x2 | 1, "TRE_Z", pd, individual = "name",
                             integration = "quadrature", keep_objective = TRUE,
                             maxit.bobyqa = 30, maxit.optim = 1, optHessian = FALSE))
  r <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pd, individual = "name",
                             integration = "quadrature", keep_objective = TRUE,
                             maxit.bobyqa = 30, maxit.optim = 1, optHessian = FALSE))
  expect_identical(z$integration, "quadrature")
  ## TRE_Z: (sigv, sigr, beta, delta) with sigma_u = sqrt(exp(delta)).
  ## TRE:   (lambda, sigma, sigr, beta).
  xz <- z$opt$par
  nb <- 3L
  sv <- xz[1]; su <- sqrt(exp(xz[length(xz)]))
  xr <- c(su / sv, sqrt(sv^2 + su^2), xz[2], xz[3:(2 + nb)])
  expect_equal(z$objective(xz, per_obs = TRUE), r$objective(xr, per_obs = TRUE),
               tolerance = 1e-12)
})

test_that("quadrature is refused outside TRE/TRE_Z, and the default is unchanged", {
  pd <- pd_tre(N = 20, t = 4, rand = 3)
  expect_error(suppressWarnings(psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name",
                                     estimator = "sml", integration = "quadrature")),
               "two-dimensional")
  expect_error(suppressWarnings(psfm(y_gtre ~ x1 + x2, "TFE", pd, individual = "name",
                                     integration = "quadrature")),
               "implemented for model_name")
  skip_on_cran()
  f <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pd, individual = "name",
                             halton_num = 20, rand.gtre = 7L, maxit.bobyqa = 20,
                             maxit.optim = 1, optHessian = FALSE))
  expect_identical(f$integration, "simulate")
})
