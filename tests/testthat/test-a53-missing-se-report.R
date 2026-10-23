## A53 item 2. A ttsfm() fit whose Hessian is indefinite loses the standard
## errors of the negative-variance coordinates and used to say nothing about
## it, after minutes of fitting. Driven on hand-built matrices rather than
## through a fit: the fits that reach this state sit on flat ridges, which is
## what does not reproduce across platforms. The live reproduction is the slow
## test at the bottom.

mk_h <- function(ev) {
  set.seed(11)
  Q <- qr.Q(qr(matrix(stats::rnorm(length(ev)^2), length(ev), length(ev))))
  Q %*% diag(ev) %*% t(Q)
}

test_that("nothing is reported when every standard error is finite", {
  expect_null(sfa:::.se_missing_report(c(a = 1, b = 2), c("a", "b"), diag(2)))
})

test_that("the missing parameters are named, with the conditioning", {
  h <- mk_h(c(30578.7, 300, 100, 10, 1, -0.00311802))
  se <- c(5.03, 1.48, 1.99, NaN, 0.51, NaN)
  nms <- c("(Intercept)", "x1", "x2", "sigma_v", "sigma_u", "sigma_w")
  r <- sfa:::.se_missing_report(se, nms, h, "ttsfm(model_name = \"TTHN\"): ")
  expect_equal(r$parameters, c("sigma_v", "sigma_w"))
  expect_equal(r$n_missing, 2L)
  expect_equal(r$n_par, 6L)
  expect_true(r$indefinite)
  expect_equal(r$min_eigenvalue, -0.00311802, tolerance = 1e-9)
  expect_equal(r$max_eigenvalue, 30578.7, tolerance = 1e-9)
  expect_match(r$message, "^ttsfm\\(model_name = \"TTHN\"\\): ")
  expect_match(r$message, "sigma_v, sigma_w", fixed = TRUE)
  expect_match(r$message, "not positive definite", fixed = TRUE)
  expect_match(r$message, "-0.00312", fixed = TRUE)
})

test_that("the message does not decide between the two readings", {
  ## A boundary optimum has no reason to have a positive definite Hessian in
  ## the full parameter space, so "not a maximum" must not be asserted alone.
  r <- sfa:::.se_missing_report(c(1, NaN), c("a", "b"), mk_h(c(10, -1)))
  expect_match(r$message, "EITHER a point that is not a maximum OR a boundary",
    fixed = TRUE
  )
  expect_match(r$message, "tt_boundary_report()", fixed = TRUE)
})

test_that("an unreadable Hessian is reported as unavailable, not as indefinite", {
  r <- sfa:::.se_missing_report(c(1, NA_real_), c("a", "b"), NULL)
  expect_true(is.na(r$indefinite))
  expect_true(is.na(r$min_eigenvalue))
  expect_match(r$message, "could not be read", fixed = TRUE)
  expect_false(grepl("not positive definite", r$message, fixed = TRUE))
  r2 <- sfa:::.se_missing_report(c(1, NA_real_), c("a", "b"),
    matrix(c(1, NA, NA, 1), 2, 2)
  )
  expect_true(is.na(r2$indefinite))
})

test_that("a positive definite but near-singular Hessian gets its own wording", {
  r <- sfa:::.se_missing_report(c(1, NaN), c("a", "b"), mk_h(c(1e12, 1e-14)))
  expect_false(r$indefinite)
  expect_match(r$message, "close enough to singular", fixed = TRUE)
  expect_false(grepl("not positive definite", r$message, fixed = TRUE))
})

test_that("unnamed standard errors fall back to positional labels", {
  r <- sfa:::.se_missing_report(c(1, NaN), character(0), diag(c(1, -1)))
  expect_equal(r$parameters, "[2]")
})

test_that("print() reports the flag, because a warning does not survive saveRDS", {
  x <- structure(
    list(
      out = matrix(1, 1, 1), formula = y ~ x, total_time = 1,
      opt = list(value = 1), model_name = "TTHN",
      se_unresolved = list(parameters = c("sigma_v", "sigma_w"), n_missing = 2L,
                           n_par = 6L, indefinite = TRUE)
    ),
    class = "sfareg"
  )
  out <- paste(capture.output(print(x)), collapse = " ")
  expect_match(out, "no standard error could be computed for sigma_v, sigma_w",
    fixed = TRUE
  )
  expect_match(out, "2 of 6 parameters", fixed = TRUE)
  expect_match(out, "not positive definite", fixed = TRUE)
})

test_that("a TTHN fit with an indefinite Hessian warns and carries the flag", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("SFA_TEST_SLOW")),
    "one TTHN fit is ~400s; set SFA_TEST_SLOW to run"
  )
  ## Seed 3002 of A53's design: the fit's own Hessian has one negative
  ## eigenvalue (-0.00312) and two of six standard errors come back NaN.
  d <- data_gen_cs(N = 300, rand = 3002, sig_u = .5, sig_v = .2, cons = 5,
                   beta1 = 1.5, beta2 = 2, a = 5, mu = .1)
  expect_warning(
    f <- ttsfm(y_tthn ~ x1 + x2, model_name = "TTHN", data = d),
    "no standard error could be computed"
  )
  expect_false(is.null(f$se_unresolved))
  expect_true(f$se_unresolved$indefinite)
  expect_equal(f$se_unresolved$n_missing, sum(!is.finite(f$out[, "st_err"])))
  expect_match(paste(capture.output(print(f)), collapse = " "),
    "no standard error could be computed", fixed = TRUE
  )
})
