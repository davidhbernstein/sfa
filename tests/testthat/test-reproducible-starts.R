## pcs_c() drew its starting values with unseeded runif(), so FD, TFE and
## TFE_WMLE fits were not reproducible run to run -- two fits of identical data
## in one session differed by up to 2.8e-4 relative on FD -- and every such fit
## silently advanced the caller's random-number stream.

.repro_panel <- function() {
  as.data.frame(data_gen_p(t = 5, N = 40, rand = 8, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
}

test_that("pcs_c() is deterministic and leaves the caller's RNG state alone", {
  y <- -abs(rnorm(200)) + rnorm(200, 0, 0.3)
  set.seed(1); before <- .Random.seed
  a <- pcs_c(y)[[1]]$par
  expect_identical(.Random.seed, before)
  set.seed(999)
  b <- pcs_c(y)[[1]]$par
  expect_identical(a, b)
})

test_that("FD fits are reproducible and do not advance the caller's RNG", {
  skip_on_cran()
  d <- .repro_panel()
  set.seed(1); before <- .Random.seed
  f1 <- suppressWarnings(psfm(y_fd ~ x_fd | z_fd, model_name = "FD", data = d,
    individual = "name", time = "year"))
  expect_identical(.Random.seed, before)
  set.seed(999)
  f2 <- suppressWarnings(psfm(y_fd ~ x_fd | z_fd, model_name = "FD", data = d,
    individual = "name", time = "year"))
  expect_identical(coef(f1), coef(f2))
})

test_that("TFE_WMLE fits are reproducible", {
  skip_on_cran()
  d <- .repro_panel()
  set.seed(1)
  f1 <- suppressWarnings(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d,
    individual = "name", time = "year"))
  set.seed(999)
  f2 <- suppressWarnings(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d,
    individual = "name", time = "year"))
  expect_identical(coef(f1), coef(f2))
})
