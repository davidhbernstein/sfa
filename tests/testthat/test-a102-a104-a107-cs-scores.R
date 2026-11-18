## A102, A104, A107(a): cross-sectional efficiency scores that disagreed with
## their own model. Each test pins arithmetic or an invariant the fix
## guarantees, not where an optimizer lands.

test_that("cols and cmle score a cost frontier on the production-oriented residual (A107a)", {
  ## -y is a cost frontier with the frontier negated and the same u and v.
  ## Both estimators fit Yc = inefdec_n * Y, so the cost fit of -y and the
  ## production fit of y see the SAME data internally: their scores must be
  ## the same. Before the fix the cost residual was flipped a second time.
  d <- cs_small(N = 300)
  dc <- d
  dc$y_pcs <- -dc$y_pcs
  for (est in c("cols", "cmle")) {
    p <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN",
      estimator = est))
    cst <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = dc, model_name = "NHN",
      estimator = est, inefdec = FALSE))
    skip_if(all(is.na(p$exp_u_hat)), "sigma_u collapsed under wrong skew")
    expect_equal(as.numeric(cst$exp_u_hat), as.numeric(p$exp_u_hat),
      tolerance = 1e-6, label = est)
  }
})

test_that("sfm('NR') stores E[u | e], checked against a posterior built by hand (A104)", {
  d <- cs_small(N = 300)
  f <- suppressWarnings(sfm(y_pcs_r ~ x1 + x2, data = d, model_name = "NR"))
  expect_false(is.null(f$u_hat))
  b <- f$out[, "par"]
  sv <- b[[1]]; su <- b[[2]]; beta <- b[-(1:2)]
  e <- as.numeric(d$y_pcs_r - cbind(1, d$x1, d$x2) %*% beta)
  ## The Rayleigh density the package uses, f(u) = 2u/su^2 exp(-u^2/su^2).
  ## exp_u_hat must come out of the same posterior, which certifies the density.
  post <- function(i, g) {
    k <- function(u) dnorm(e[i] + u, 0, sv) * 2 * u / su^2 * exp(-u^2 / su^2)
    integrate(function(u) g(u) * k(u), 0, Inf, rel.tol = 1e-11)$value /
      integrate(k, 0, Inf, rel.tol = 1e-11)$value
  }
  idx <- c(order(e)[1:3], order(-e)[1:3], 150L)
  ## The test points include both tails.
  want_u <- vapply(idx, post, numeric(1), g = identity)
  want_b <- vapply(idx, post, numeric(1), g = function(u) exp(-u))
  expect_equal(as.numeric(f$exp_u_hat[idx]), want_b, tolerance = 1e-7)
  expect_equal(as.numeric(f$u_hat[idx]), want_u, tolerance = 1e-7)
  ## And it is not -log(exp_u_hat), which is what sfma() used to average.
  expect_gt(max(abs(f$u_hat + log(f$exp_u_hat))), 1e-3)
})

test_that("ivsfm() and selsfm() report E[exp(-u) | e] under type = 'bc' (A102)", {
  set.seed(1)
  n <- 600
  w1 <- rnorm(n); w2 <- rnorm(n); x1 <- rnorm(n); eta <- rnorm(n)
  v <- 0.5 * (0.6 * eta + 0.8 * rnorm(n))
  x2 <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
  d <- data.frame(y = 0.5 + 0.8 * x1 - 0.6 * x2 + v - abs(rnorm(n)),
    x1 = x1, x2 = x2, w1 = w1, w2 = w2)
  for (mn in c("IVLIML", "IVCF", "C2SLS")) {
    f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2,
      instruments = ~ w1 + w2, data = d, model_name = mn))
    bc <- as.numeric(efficiency(f, type = "bc"))
    jl <- as.numeric(efficiency(f, type = "jlms"))
    ## Jensen: E[exp(-u)] >= exp(-E[u]), strictly for a non-degenerate
    ## posterior. "bc" used to return exp(-jlms) itself.
    expect_true(all(bc >= jl - 1e-12), label = mn)
    expect_gt(max(bc - jl), 1e-3, label = mn)
  }

  ## selsfm("greene"): over the same simulated weights, so Jensen holds draw
  ## by draw whatever Nsim is.
  set.seed(2)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n)
  sel <- as.numeric(0.2 + 0.7 * z1 - 0.5 * z2 + w > 0)
  ys <- 0.5 + 0.8 * x1 - 0.3 * x2 + 0.5 * (0.6 * w + 0.8 * rnorm(n)) - abs(rnorm(n))
  ys[sel == 0] <- NA
  ds <- data.frame(y = ys, x1 = x1, x2 = x2, z1 = z1, z2 = z2, d = sel)
  g <- suppressWarnings(selsfm(d ~ z1 + z2, y ~ x1 + x2, ds, Nsim = 60, seed = 1))
  bc <- as.numeric(efficiency(g, type = "bc"))
  jl <- as.numeric(efficiency(g, type = "jlms"))
  expect_true(all(bc >= jl - 1e-12))
  expect_gt(max(bc - jl), 1e-3)
})
