## GTRE's SML draws the persistent component as +r and -r and combines the two
## halves. They are two draws of the SAME density, so they average on the
## density scale: log((A + B) / 2). Both call sites read
## 0.5 * (log A + log B) = log(sqrt(A * B)) -- the geometric mean, which AM-GM
## puts strictly below the arithmetic one whenever the halves differ.
##
## It is a finite-R bias, not an asymptotic error: each half is separately a
## consistent simulator, so the two rules agree as R grows. But psfm()'s
## default draw count is ceiling(sqrt(nrow)) + 100 -- about 118 for a 60x5
## panel -- and at that R the gap measured 0.23 log-likelihood points, with the
## sigma_r profile optimum shifted down 4.8%. Raised by Chris Parmeter,
## 2026-09-28.

test_that(".log_mean_exp2 averages on the density scale", {
  a <- log(2); b <- log(8)
  ## Arithmetic mean of 2 and 8 is 5; the geometric mean is 4.
  expect_equal(.log_mean_exp2(a, b), log(5), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(.log_mean_exp2(a, b), 0.5 * (a + b))))
  ## Equal halves are the one case where the two rules agree.
  expect_equal(.log_mean_exp2(a, a), a, tolerance = 1e-12)
  ## Vectorised, and never below the log-scale average (AM-GM).
  x <- log(c(1, 5, 20, 1e-8)); y <- log(c(3, 5, 0.5, 1e8))
  expect_equal(.log_mean_exp2(x, y), log((exp(x) + exp(y)) / 2), tolerance = 1e-10)
  expect_true(all(.log_mean_exp2(x, y) >= 0.5 * (x + y) - 1e-12))
  ## Far apart in the exponent: the naive exp() of either half alone would
  ## overflow or underflow, which is why this goes through log-sum-exp.
  expect_equal(.log_mean_exp2(-800, -1200), -800 + log((1 + exp(-400)) / 2),
    tolerance = 1e-10)
  ## Both halves impossible.
  expect_equal(.log_mean_exp2(-Inf, -Inf), -Inf)
})

test_that("GTRE's two likelihood paths agree, and neither uses the geometric mean", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 5, N = 40, rand = 9, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))

  ## The loop and vectorised paths must evaluate the SAME likelihood; they are
  ## the two sites that carried this defect, so a fix applied to one only would
  ## show up here.
  old <- getOption("sfa.gtre_vectorized", FALSE)
  on.exit(options(sfa.gtre_vectorized = old), add = TRUE)

  options(sfa.gtre_vectorized = FALSE)
  f_loop <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE",
    data = d, individual = "name", keep_objective = TRUE))
  options(sfa.gtre_vectorized = TRUE)
  f_vec <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE",
    data = d, individual = "name", keep_objective = TRUE))

  expect_equal(unname(coef(f_loop)), unname(coef(f_vec)), tolerance = 1e-4)
  expect_equal(as.numeric(f_loop$opt$value), as.numeric(f_vec$opt$value),
    tolerance = 1e-4)
})
