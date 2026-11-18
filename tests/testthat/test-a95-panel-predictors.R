## A95, A98, A99, A107: panel efficiency scores that disagreed with their own
## model. Each test pins arithmetic or an invariant the fix guarantees, not
## where an optimizer lands.

test_that("GTRE (fiml) scores the UNCLAMPED residual (A95)", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 40)
  f <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = d, model_name = "GTRE",
    individual = "name", time = "year", estimator = "fiml"))
  b <- coef(f)
  np <- length(b)
  dd <- f$data
  eps <- dd$y_gtre - b[[1]] - b[["x1"]] * dd$x1 - b[["x2"]] * dd$x2
  ## Some residuals are positive, so a clamp would have changed something.
  expect_gt(mean(eps > 0), 0.01)
  e_i <- split(eps, factor(dd$name, levels = unique(dd$name)))
  te <- .gtre_csn_te(e_i, rep(6L, 40L), sig_u = b[[np]], sig_v = b[[np - 2]],
    sig_r = b[[np - 3]], sig_h = b[[np - 1]])
  expect_equal(as.numeric(f$U), as.numeric(te$U), tolerance = 1e-10)
  expect_equal(as.numeric(f$H), as.numeric(te$H), tolerance = 1e-10)
  ## The clamped residual gives a different answer on the same parameters.
  te_clamped <- .gtre_csn_te(lapply(e_i, function(z) pmin(z, 0)), rep(6L, 40L),
    b[[np]], b[[np - 2]], b[[np - 3]], b[[np - 1]])
  expect_gt(max(abs(te_clamped$U - te$U)), 1e-3)
})

test_that("GTRE (fiml) on a cost frontier mirrors the production fit (A107)", {
  skip_on_cran()
  ## -y is a cost frontier with every frontier coefficient negated and the
  ## same scales, so the cost fit must report exactly that -- including the
  ## intercept, which was reported with the wrong sign.
  d <- panel_small(t = 6, N = 40)
  dc <- d
  dc$y_gtre <- -dc$y_gtre
  p <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = d, model_name = "GTRE",
    individual = "name", time = "year", estimator = "fiml"))
  cst <- suppressWarnings(psfm(y_gtre ~ x1 + x2, data = dc, model_name = "GTRE",
    individual = "name", time = "year", estimator = "fiml", inefdec = FALSE))
  fr <- c("(Intercept)", "x1", "x2")
  sc <- c("sigr", "sigv", "sigh", "sigu")
  expect_equal(unname(coef(cst)[fr]), unname(-coef(p)[fr]), tolerance = 1e-3)
  expect_equal(unname(coef(cst)[sc]), unname(coef(p)[sc]), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(cst)), as.numeric(logLik(p)), tolerance = 1e-6)
  ## And the scores are of the same posterior.
  expect_equal(as.numeric(cst$U), as.numeric(p$U), tolerance = 1e-3)
  expect_equal(as.numeric(cst$H), as.numeric(p$H), tolerance = 1e-3)
})

test_that("FD's exp_u_hat is E[exp(-u) | e], not exp(-E[u | e]) (A99)", {
  skip_on_cran()
  ## A non-null design (mu = 0.5, delta = 0.3), where the posterior of u*_i
  ## is not degenerate and the two quantities differ visibly.
  d <- as.data.frame(data_gen_p(t = 6, N = 40, rand = 7, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5,
    mu = 0.5, delta_fd = 0.3, fd_draw = "truncated"))
  fd <- suppressWarnings(psfm(y_fd ~ x_fd | z_fd, model_name = "FD",
    data = d, individual = "name", time = "year"))
  ## Jensen: E[exp(-u)] >= exp(-E[u]) for every observation, strictly unless
  ## the posterior is degenerate. The old slot WAS exp(-u_hat).
  expect_true(all(fd$exp_u_hat >= exp(-fd$u_hat) - 1e-12))
  expect_gt(max(fd$exp_u_hat - exp(-fd$u_hat)), 1e-3)
  expect_true(all(fd$exp_u_hat > 0 & fd$exp_u_hat <= 1))
})

test_that("TFE_WMLE's exp_u_hat is not tied at a clamp (A98)", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 40)
  w <- suppressWarnings(psfm(y_tfe ~ x1 + x2, model_name = "TFE_WMLE",
    data = d, individual = "name", time = "year"))
  ## Under the clamp every observation with a positive residual tied at one
  ## value, the maximum. Battese-Coelli is strictly monotone in the residual.
  expect_lte(sum(abs(w$exp_u_hat - max(w$exp_u_hat)) < 1e-12), 1L)
})
