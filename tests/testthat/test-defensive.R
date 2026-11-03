## Package hygiene: RNG state, and argument validation on the exported fitters.

test_that("seeded functions restore the caller's RNG stream", {
  ## A package function that calls set.seed() overwrites .Random.seed in the
  ## global environment, so the user's own seed stops controlling their session
  ## the moment they fit a model. Demonstrated before the fix: with the caller
  ## holding set.seed(42), rnorm(3) gave 1.37096 -0.56470 0.36313 alone but
  ## 0.99032 0.11227 1.14963 after an intervening data_gen_cs(rand = 7).
  args <- list(N = 10, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
               beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)

  set.seed(42); ref <- stats::rnorm(3)
  set.seed(42); invisible(do.call(data_gen_cs, args)); got <- stats::rnorm(3)
  expect_equal(got, ref)

  set.seed(42); ref_p <- stats::rnorm(3)
  set.seed(42)
  invisible(data_gen_p(t = 3, N = 10, rand = 7, sig_u = 1, sig_v = 0.3,
                       sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                       beta1 = 0.5, beta2 = 0.5))
  expect_equal(stats::rnorm(3), ref_p)
})

test_that("an unseeded call still advances the stream, as ordinary code does", {
  ## Restoring must be conditional on our having hijacked the seed. With no
  ## seed supplied the draws should come from, and advance, the caller's stream.
  set.seed(11)
  a <- data_gen_cs(N = 8, rand = NULL, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)$y_pcs
  b <- data_gen_cs(N = 8, rand = NULL, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)$y_pcs
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("no .Random.seed is left behind when none existed", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(data_gen_cs(N = 8, rand = 3, sig_u = 1, sig_v = 0.3, cons = 0.5,
                        beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("seeding still makes the generators reproducible", {
  ## The restore must not change what is produced, only what is left behind.
  args <- list(N = 50, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
               beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_equal(do.call(data_gen_cs, args), do.call(data_gen_cs, args))
})

test_that("exported fitters reject malformed arguments with a named message", {
  d <- cs_small(N = 60)
  for (fn in c("sfm", "zsfm", "ttsfm")) {
    f <- get(fn, envir = asNamespace("sfa"))
    expect_error(f("y ~ x", data = d), "`formula` must be a formula")
    expect_error(f(y_pcs ~ x1 + x2), "`data` must be supplied")
    expect_error(f(y_pcs ~ x1 + x2, data = as.matrix(d)), "`data` must be a data.frame")
    expect_error(f(y_pcs ~ x1 + x2, data = d[0, ]), "at least 1 is required")
  }
  expect_error(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN", maxit.optim = -1),
               "`maxit.optim` must be a single positive number")
  expect_error(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN", optHessian = "yes"),
               "`optHessian` must be TRUE or FALSE")
  expect_error(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN", inefdec = NA),
               "`inefdec` must be TRUE or FALSE")
})

test_that("every fitter refuses an orientation that is not TRUE or FALSE (A70)", {
  ## isTRUE() maps every non-TRUE value to FALSE, so before this check psfm()
  ## and esfm() fitted a COST frontier for inefdec = "yes" and recorded it on
  ## the fit as though it had been asked for.
  d <- cs_small(N = 60)
  p <- data.frame(id = rep(1:12, each = 5), x1 = d$x1, y = d$y_pcs)
  fits <- list(
    sfm   = function(v) sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN", inefdec = v),
    zsfm  = function(v) zsfm(y_pcs ~ x1 + x2, data = d, inefdec = v),
    lcsfm = function(v) lcsfm(y_pcs ~ x1 + x2, data = d, inefdec = v),
    ttsfm = function(v) ttsfm(y_pcs ~ x1 + x2, data = d, inefdec = v),
    psfm  = function(v) psfm(y ~ x1, data = p, individual = "id",
                             model_name = "TRE", inefdec = v),
    esfm  = function(v) esfm(y_pcs ~ x1 + x2, data = d, inefdec = v)
  )
  for (fn in names(fits)) {
    for (v in list(NA, "yes", 1, c(TRUE, TRUE))) {
      expect_error(fits[[fn]](v), "`inefdec` must be TRUE or FALSE",
                   info = paste(fn, deparse(v)))
    }
  }
})

test_that("psfm() validates its arguments as the other fitters do (A70)", {
  d <- cs_small(N = 60)
  p <- data.frame(id = rep(1:12, each = 5), x1 = d$x1, y = d$y_pcs)
  expect_error(psfm("y ~ x1", data = p, individual = "id"),
               "`formula` must be a formula")
  expect_error(psfm(y ~ x1, individual = "id"), "`data` must be supplied")
  expect_error(psfm(y ~ x1, data = as.matrix(p), individual = "id"),
               "`data` must be a data.frame")
  expect_error(psfm(y ~ x1, data = p, individual = "id", maxit.optim = -1),
               "`maxit.optim` must be a single positive number")
  expect_error(psfm(y ~ x1, data = p, individual = "id", PSopt = "no"),
               "`PSopt` must be TRUE or FALSE")
})

test_that("validation accepts tibble-like and pdata.frame-like inputs", {
  ## inherits(data, "data.frame") rather than an identity check on class(),
  ## so subclasses keep working.
  d <- cs_small(N = 60)
  class(d) <- c("tbl_df", "tbl", "data.frame")
  expect_error(sfa:::.validate_sfa_call(y_pcs ~ x1 + x2, d, "sfm"), NA)
})

test_that("a singular Hessian costs only the standard errors, not the fit", {
  ## sfm() used to fill the out matrix with try(..., silent = T), which both
  ## used the rebindable T and could leave a row half-written. The contract is
  ## that point estimates always survive and the SE/t rows degrade to NA.
  skip_on_cran()
  d   <- cs_small(N = 80)
  fit <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN"))
  expect_true(all(is.finite(fit$out[, "par"])))
  expect_equal(nrow(fit$out), length(coef(fit)))
  expect_true(is.numeric(fit$out[, "st_err"]))
})
