## A62. optim() can return convergence 0 at a point that is not a maximum, and
## sfm() then reports standard errors for it. .stationarity_polish() checks the
## Newton decrement at the returned point and polishes when it predicts a gain.
## Driven on synthetic objectives rather than through a fit: the fits that reach
## this state sit on flat ridges, which is exactly what does not reproduce
## across platforms (the lesson from #54, #57 and #59).

mk_opt <- function(fn, par, hess) {
  list(par = par, value = fn(par), convergence = 0L, message = NULL,
       hessian = hess)
}

test_that("a point at the minimum is left alone", {
  fn <- function(x) sum((x - c(1, 2))^2)
  o <- mk_opt(fn, c(1, 2), diag(2, 2))
  r <- sfa:::.stationarity_polish(fn, o, lower = c(-Inf, -Inf))
  expect_false(r$info$polished)
  expect_lt(abs(r$info$decrement), 1e-8)
  expect_identical(r$opt, o)
})

test_that("a point off the minimum is polished to it", {
  fn <- function(x) sum((x - c(1, 2))^2)
  o <- mk_opt(fn, c(3, -1), diag(2, 2))
  r <- sfa:::.stationarity_polish(fn, o, lower = c(-Inf, -Inf))
  ## 0.5 g'H^-1 g for this quadratic is exactly the gap, 4 + 9 = 13.
  expect_equal(r$info$decrement, 13, tolerance = 1e-6)
  expect_true(r$info$polished)
  expect_equal(r$opt$par, c(1, 2), tolerance = 1e-5)
  expect_equal(r$info$gain, 13, tolerance = 1e-6)
  expect_lt(abs(r$info$decrement_after), 1e-4)
  expect_match(r$opt$message, "polished by nlminb")
})

test_that("a ridge to a positive parameter's infinity is walked out along", {
  ## The supremum is x1 -> Inf, as for NTN's lambda on a sigma_v collapse.
  fn <- function(x) exp(-x[1]) + (x[2] - 1)^2
  o <- mk_opt(fn, c(2, 1), diag(c(exp(-2), 2)))
  r <- sfa:::.stationarity_polish(fn, o, lower = c(.Machine$double.eps, -Inf))
  expect_true(r$info$polished)
  expect_gt(r$opt$par[1], 10)
  expect_lt(r$opt$value, 1e-4)
  expect_equal(r$opt$par[2], 1, tolerance = 1e-4)
})

test_that("a singular Hessian is not acted on", {
  fn <- function(x) sum((x - c(1, 2))^2)
  o <- mk_opt(fn, c(3, -1), matrix(0, 2, 2))
  r <- sfa:::.stationarity_polish(fn, o, lower = c(-Inf, -Inf))
  expect_true(is.na(r$info$decrement))
  expect_false(r$info$polished)
  expect_identical(r$opt, o)
})

test_that("a polish that cannot improve leaves the fit unchanged", {
  ## The decrement is computed from the REPORTED Hessian, which can be wrong;
  ## a wrong one must not move a point that is in fact the minimum.
  fn <- function(x) sum((x - c(1, 2))^2)
  o <- mk_opt(fn, c(1, 2), diag(2, 2))
  fn_g <- fn
  o$hessian <- diag(1e-12, 2)
  r <- sfa:::.stationarity_polish(fn_g, o, lower = c(-Inf, -Inf))
  expect_false(r$info$polished)
  expect_identical(r$opt$par, c(1, 2))
})

test_that("an interior sfm() fit carries $stationarity and is not polished", {
  d <- data_gen_cs(300, 1, sig_u = 1, sig_v = 0.3, cons = 0.5, beta1 = 0.5,
    beta2 = 0.5, a = 5, mu = 0.5)
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN"))
  expect_type(f$stationarity, "list")
  expect_false(f$stationarity$polished)
  ## Loose on purpose: a decrement this small is many orders below the 1e-4
  ## threshold on any platform, so the assertion does not pin a ridge.
  expect_lt(abs(f$stationarity$decrement), 1e-6)
})
