## log D_nu(z), the parabolic cylinder function behind NG and NNAK. Computed in
## pure R since gap A25 retired gsl; see notes/code_history/matrix_utils.md.

test_that("log D_-1 matches its closed form wherever the value is exact", {
  ## D_-1(z) = exp(z^2/4) sqrt(pi/2) erfc(z/sqrt(2)): an exact reference needing
  ## no special-function package, across the series, subtraction and peak
  ## branches. Below z = -sqrt(1400) the old clipped value is reproduced instead
  ## (next tests).
  z <- c(-37.4, -37, -12, -10, -5, -1, 0, 0.3, 0.5, 0.51, 1, 3, 6, 12, 20, 60, 300)
  exact <- z^2 / 4 + 0.5 * log(pi / 2) + log(2) + pnorm(-z, log.p = TRUE)
  expect_equal(.log_pcf(-1, z), exact, tolerance = 1e-10)
})

test_that(".log_pcf agrees with gsl's hypergeometric route where that works", {
  skip_if_not_installed("gsl")
  ## The implementation it replaced: U for z > 0.5, a difference of two 1F1
  ## terms at or below it. Compared only where those neither overflow nor cancel.
  worst <- 0
  for (nu in c(-0.2, -0.5, -1, -2, -4.11, -8, -12)) {
    for (z in c(-8, -3, -1, 0, 0.3, 0.6, 1, 2, 5, 10, 20)) {
      a <- if (z > 0.5) {
        u <- tryCatch(gsl::hyperg_U(-nu / 2, 0.5, z^2 / 2), error = function(e) NA_real_)
        (nu / 2) * log(2) - z^2 / 4 + log(u)
      } else {
        q <- z^2 / 2
        br <- gsl::hyperg_1F1(-nu / 2, 0.5, q) / gamma((1 - nu) / 2) -
          sqrt(2) * z * gsl::hyperg_1F1((1 - nu) / 2, 1.5, q) / gamma(-nu / 2)
        if (!is.finite(br) || br <= 0) NA_real_ else (nu / 2) * log(2) + 0.5 * log(pi) - z^2 / 4 + log(br)
      }
      if (!is.finite(a)) next
      worst <- max(worst, abs(.log_pcf(nu, z) - a))
    }
  }
  expect_lt(worst, 1e-8)
})

test_that("it is finite above the clip, and never NaN anywhere", {
  G <- expand.grid(nu = -c(0.01, 0.05, 0.5, 1, 4, 16, 64), z = c(-37.4, -10, -1, 0, 0.5, 2, 12, 40, 300))
  expect_silent(v <- .log_pcf(G$nu, G$z))
  expect_true(all(is.finite(v)))
  ## Below the clip: finite, or +Inf where the old code overflowed -- never NaN.
  w <- .log_pcf(-c(0.5, 4, 16, 64), c(-40, -40, -300, -300))
  expect_false(anyNA(w))
  expect_true(all(is.finite(w) | w == Inf))
})

test_that("below the clip it reproduces the gsl code it replaced, overflow included", {
  skip_if_not_installed("gsl")
  ## NG's optimum can sit against this edge, so the old (clipped, understated)
  ## value is reproduced rather than corrected; see notes/code_history.
  old <- function(p, z) {
    q <- pmin(z^2 / 2, 700)
    br <- suppressWarnings(gsl::hyperg_1F1(p / 2, 0.5, q) / gamma((1 + p) / 2) -
      sqrt(2) * z * gsl::hyperg_1F1((1 + p) / 2, 1.5, q) / gamma(p / 2))
    -(p / 2) * log(2) + 0.5 * log(pi) - z^2 / 4 + log(br)
  }
  G <- expand.grid(p = c(0.02, 0.5, 1.3, 2.05, 4.11, 10.5), z = c(-38, -45, -60, -100, -300))
  o <- old(G$p, G$z)
  n <- .log_pcf(-G$p, G$z)
  expect_identical(is.finite(n), is.finite(o))
  expect_equal(n[is.finite(o)], o[is.finite(o)], tolerance = 1e-10)
})

test_that("a large shape returns a finite, decreasing value", {
  for (nu in c(-32, -64)) {
    v <- .log_pcf(nu, c(1, 5, 20, 40))
    expect_true(all(is.finite(v)), info = paste("nu =", nu))
    expect_true(all(diff(v) < 0), info = paste("nu =", nu))
  }
})

test_that("the integrate() reference still works for small shapes", {
  expect_silent(v <- .log_pcf_integral(-0.5, c(0.6, 1, 5)))
  expect_true(all(is.finite(v)))
  expect_equal(.log_pcf(-0.5, c(0.6, 1, 5)), v, tolerance = 1e-8)
})

test_that("NNAK fits at shapes that used to kill it", {
  skip_on_cran()
  for (m in c(1, 4)) {
    d <- data_gen_cs(N = 400, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
      beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5, m_nak = m
    )
    f <- try(suppressWarnings(
      sfm(y_pcs_nak ~ x1 + x2, model_name = "NNAK", data = d)
    ), silent = TRUE)
    expect_false(inherits(f, "try-error"), info = paste("m_nak =", m))
    expect_true(is.finite(f$opt$value), info = paste("m_nak =", m))
  }
})

test_that("log D_nu accepts a VECTOR order, one per observation", {
  ## This is what lets the gamma/Nakagami shape depend on covariates (G5).
  z <- c(-2, -0.3, 0.6, 2, 8)
  expect_equal(.log_pcf(-4, z), .log_pcf(rep(-4, length(z)), z))
  nu <- c(-1, -2, -4, -8, -16)
  v <- .log_pcf(nu, z)
  ref <- vapply(seq_along(z), function(i) .log_pcf(nu[i], z[i]), numeric(1))
  expect_equal(v, ref)
  expect_true(all(is.finite(v)))
})

test_that("an NG fit is never worse than the start its optimizer was handed", {
  skip_on_cran()
  ## On this sample the three optimizer stages took NG's polished multistart
  ## point (log-likelihood -259.6) and returned -449.5 (gap A25).
  set.seed(4)
  x1 <- rnorm(300); x2 <- rnorm(300)
  y <- 2 + 0.5 * x1 + 0.5 * x2 + rnorm(300, 0, 0.3) - rgamma(300, shape = 2, scale = 0.4)
  f <- suppressWarnings(sfm(y ~ x1 + x2, model_name = "NG", data = data.frame(y, x1, x2)))
  expect_gte(-f$opt$value, f$ng_starts$best - 1e-6)
  expect_gt(-f$opt$value, -300)
})

test_that("NG's fitted values are unchanged", {
  skip_on_cran()
  ## A REGRESSION PIN. Removing the old accidental barrier moved NG from
  ## (sigv 0.185, x1 0.560) to (sigv 0.000, x1 0.775) against a true x1 of 0.5
  ## -- at a HIGHER likelihood, so no likelihood check would catch it. The old
  ## values below the clip are reproduced for that reason; these must not move.
  d <- data_gen_cs(N = 400, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5
  )
  f <- suppressWarnings(sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d))
  p <- f$out[, "par"]

  expect_equal(unname(p[["sigv"]]), 0.185, tolerance = 0.02)
  expect_equal(unname(p[["sigu"]]), 0.373, tolerance = 0.02)
  expect_equal(unname(p[["x1"]]), 0.560, tolerance = 0.02)
  expect_gt(p[["sigv"]], 0.05)
})
