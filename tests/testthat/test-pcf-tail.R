## Issue #30. Two defects in NG's and NNAK's likelihood for z far below zero:
## .log_pcf() reproduced the gsl code's CLIPPED value below z = -37.42 (wrong
## by hundreds of log units), and the densities formed z^2/2 - eps^2/(2 sig_v^2)
## as two terms of ~1e14 each once sigma_v neared its floor. Every check below
## is against a closed form or a numerical integral, not another fitted value.

## nu = -1: D_{-1}(z) = sqrt(2 pi) exp(z^2/4) Phi(-z).
d1 <- function(z) z^2 / 4 + 0.5 * log(2 * pi) + pnorm(-z, log.p = TRUE)

test_that(".log_pcf() is exact below the old cutoff", {
  z <- c(-37.5, -42.9, -60, -100, -300)
  expect_equal(.log_pcf(-1, z), d1(z), tolerance = 1e-12)
})

test_that(".log_pcf_w() is exact for z far below zero", {
  z <- -c(0.5, 5, 29.9, 30.1, 42.9, 1e3, 1e5, 1e7)
  ## Written directly: d1(z) - z^2/4 would itself cancel at z = -1e7.
  expect_equal(.log_pcf_w(-1, z), 0.5 * log(2 * pi) + pnorm(-z, log.p = TRUE), tolerance = 1e-9)
  ## Any order, against the integral itself, on both sides of the z = -30 switch.
  ref <- function(p, z) {
    t0 <- (-z + sqrt(z^2 + 4 * p)) / 2
    c0 <- (p - 1) * log(t0) - (t0 + z)^2 / 2
    f <- function(t) exp((p - 1) * log(t) - (t + z)^2 / 2 - c0)
    v <- integrate(f, 0, t0, rel.tol = 1e-12, abs.tol = 0)$value +
      integrate(f, t0, Inf, rel.tol = 1e-12, abs.tol = 0)$value
    -lgamma(p) + log(v) + c0
  }
  for (p in c(0.3, 2.3, 20)) {
    for (z in c(-29, -31, -200)) {
      expect_equal(.log_pcf_w(-p, z), ref(p, z), tolerance = 1e-10,
        info = paste("p =", p, "z =", z))
    }
  }
})

test_that("NNAK at shape 1/2 is NHN, outliers included", {
  skip_on_cran()
  ## Observation 26 has z of about -43 at NHN's optimum; with the old cutoff
  ## NNAK scored -624.8 there against NHN's -404.9.
  d <- as.data.frame(data_gen_cs(N = 200, rand = 1002412, cons = 0.5, beta1 = 0.5,
    beta2 = 0.5, sig_u = 1, sig_v = 1, mu = 0.5, a = 5))
  k <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NNAK",
    keep_objective = TRUE))
  h <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, data = d, model_name = "NHN"))
  hp <- h$out[, 1]
  lam <- hp[[1]]
  sig <- hp[[2]]
  q <- unname(c(sig / sqrt(1 + lam^2), lam * sig / sqrt(1 + lam^2), 0.5, hp[3:5]))
  expect_equal(-k$objective(q), as.numeric(logLik(h)), tolerance = 1e-8)
  ## So the fit, which nests NHN, cannot end below it (the issue's symptom).
  expect_gte(as.numeric(logLik(k)), as.numeric(logLik(h)) - 1e-6)
})

test_that("NG and NNAK equal the convolution integral as sigma_v -> 0", {
  skip_on_cran()
  set.seed(30)
  n <- 40
  x1 <- rnorm(n)
  y <- 1 + 0.5 * x1 + rnorm(n, 0, 0.3) - rgamma(n, shape = 1.5, scale = 0.5)
  dd <- data.frame(y, x1)
  X <- cbind(1, x1)
  for (mod in c("NG", "NNAK")) {
    k <- suppressWarnings(sfm(y ~ x1, data = dd, model_name = mod, keep_objective = TRUE))
    lg <- if (mod == "NG") {
      function(u, m, su) dgamma(u, shape = m, scale = su, log = TRUE)
    } else {
      function(u, m, su) log(2) + m * log(m) - lgamma(m) - 2 * m * log(su) +
        (2 * m - 1) * log(u) - m * u^2 / su^2
    }
    for (sv in c(0.05, 1e-4, 1e-7)) {
      ## Frontier above every observation, so the density is not negligible.
      b <- c(max(y - 0.5 * x1) + 0.3, 0.5)
      p <- c(sv, 0.6, 1.5, b)
      e <- as.numeric(y - X %*% b)
      ref <- sum(vapply(e, function(ei) {
        f <- function(u) exp(dnorm(ei + u, 0, sv, log = TRUE) + lg(u, 1.5, 0.6))
        c0 <- -ei
        pts <- sort(unique(c(0, max(0, c0 - 40 * sv), c0, c0 + 40 * sv)))
        v <- sum(vapply(seq_len(length(pts) - 1), function(j) {
          integrate(f, pts[j], pts[j + 1], rel.tol = 1e-11, abs.tol = 0)$value
        }, numeric(1)))
        log(v + integrate(f, max(pts), Inf, rel.tol = 1e-11, abs.tol = 0)$value)
      }, numeric(1)))
      expect_equal(-k$objective(p), ref, tolerance = 1e-8,
        info = paste(mod, "sigma_v =", sv))
    }
  }
})
