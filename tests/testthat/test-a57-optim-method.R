## Gap A57. optim() applies lower/upper for "L-BFGS-B" only, and every
## opt.optim() call in this package passes bounds that hold the scale parameters
## positive. Before this guard, Method = "BFGS"/"CG"/"Nelder-Mead"/"SANN" made
## optim() warn ("bounds can only be used with method L-BFGS-B") and then search
## unbounded -- invisible on an interior optimum, and free to take a standard
## deviation to or past zero on a boundary fit.

test_that(".check_optim_method() accepts only methods that honour bounds", {
  expect_identical(.check_optim_method("L-BFGS-B"), "L-BFGS-B")
  for (m in c("BFGS", "CG", "Nelder-Mead", "SANN", "Brent")) {
    expect_error(.check_optim_method(m), "cannot be used", fixed = TRUE)
  }
  ## Shape errors are separate from the value error, and must not pass silently.
  expect_error(.check_optim_method(NULL), "single character string")
  expect_error(.check_optim_method(NA_character_), "single character string")
  expect_error(.check_optim_method(c("L-BFGS-B", "BFGS")), "single character string")
  expect_error(.check_optim_method(1L), "single character string")
})

test_that("opt.optim() refuses a method that cannot use the bounds it is given", {
  ## The invariant lives here: this function always hands optim() bounds.
  fn <- function(p, ...) sum((p - 1)^2)
  expect_error(
    opt.optim(fn, c(2, 2), lower.optim = rep(1e-7, 2), upper.optim = rep(10, 2),
              maxit.optim = 10, opt.TF = TRUE, method = "BFGS",
              optHessian = FALSE, trace = 0, verbose = FALSE),
    "cannot be used", fixed = TRUE
  )
  expect_silent(
    res <- opt.optim(fn, c(2, 2), lower.optim = rep(1e-7, 2), upper.optim = rep(10, 2),
                     maxit.optim = 10, opt.TF = TRUE, method = "L-BFGS-B",
                     optHessian = FALSE, trace = 0, verbose = FALSE)
  )
  expect_lt(res$start_feval, 2)
})

test_that("every entry point rejects the bad Method before fitting anything", {
  skip_on_cran()
  d <- data_gen_cs(N = 60, rand = 11, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  ## The error must arrive from the argument check, not from deep inside a fit,
  ## so these calls are cheap by construction.
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, Method = "BFGS"),
               "cannot be used", fixed = TRUE)
  expect_error(zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d, Method = "BFGS"),
               "cannot be used", fixed = TRUE)
  expect_error(lcsfm(y_lcm ~ x1 + x2, model_name = "LCM", data = d, Method = "BFGS"),
               "cannot be used", fixed = TRUE)
  expect_error(ttsfm(y_tthn ~ x1 + x2, model_name = "TTHN", data = d, Method = "BFGS"),
               "cannot be used", fixed = TRUE)
})

test_that("the default still fits", {
  skip_on_cran()
  d <- data_gen_cs(N = 200, rand = 11, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_s3_class(f, "sfareg")
  expect_true(all(is.finite(f$out[, "par"])))
})
