## Gap A67. Orientation was recoverable only by re-evaluating the call's
## `inefdec` expression, so an ordinary wrapper -- formula built at top level,
## fit made inside a helper with a local flag -- produced a fit on which
## efficiency(logDepVar = FALSE) and the composed CDF refused to work.

.a67_data <- function() {
  set.seed(23)
  n <- 150
  x1 <- stats::runif(n, 1, 3)
  x2 <- stats::runif(n, 1, 3)
  v <- stats::rnorm(n, 0, 0.3)
  u <- abs(stats::rnorm(n, 0, 0.6))
  data.frame(x1 = x1, x2 = x2,
    y_prod = 4 + 1.5 * x1 + 0.8 * x2 + v - u,
    y_cost = 4 + 1.5 * x1 + 0.8 * x2 + v + u
  )
}

## The formula is built HERE, so its environment is this file's, not the
## helper's -- which is why `flag` cannot be found from it later.
.a67_cost_formula <- y_cost ~ x1 + x2

test_that("a fit records the orientation it used, as a scalar logical", {
  d <- .a67_data()
  fp <- suppressWarnings(sfm(y_prod ~ x1 + x2, model_name = "NHN", data = d))
  fc <- suppressWarnings(
    sfm(y_cost ~ x1 + x2, model_name = "NHN", data = d, inefdec = FALSE)
  )
  for (f in list(fp, fc)) {
    expect_true(is.logical(f$inefdec))
    expect_length(f$inefdec, 1L)
    expect_false(is.na(f$inefdec))
  }
  expect_true(fp$inefdec)
  expect_false(fc$inefdec)
  ## The stored field and the resolver must agree in BOTH directions, or one of
  ## the two is inverted.
  expect_false(.is_cost_fit(fp))
  expect_true(.is_cost_fit(fc))
})

test_that("an inefdec the formula's environment cannot see still resolves", {
  d <- .a67_data()
  make <- function(dat) {
    flag <- FALSE
    suppressWarnings(
      sfm(.a67_cost_formula, model_name = "NHN", data = dat, inefdec = flag)
    )
  }
  f <- make(d)
  ## The premise: the call records a symbol, and evaluating it where the
  ## formula lives fails. Without that this test proves nothing.
  expect_true(is.name(f$call$inefdec))
  expect_false(exists("flag", envir = environment(f$formula), inherits = FALSE))

  expect_true(.is_cost_fit(f))
  expect_false(.sfa_inefdec(f))

  ## And the orientation reaches the quantity a user reads: a cost frontier's
  ## level-scale efficiency is f / (f + u), not the production form 1 - u / f.
  te <- efficiency(f, logDepVar = FALSE, newdata = d)
  b <- f$coefficients
  mm <- stats::model.matrix(~ x1 + x2, d)
  nm <- intersect(colnames(mm), names(b))
  fr <- as.numeric(mm[, nm, drop = FALSE] %*% b[nm])
  uu <- -log(pmax(as.numeric(f$exp_u_hat), .SFA_CONSTANTS$MIN_POSITIVE))
  expect_equal(te, fr / (fr + uu), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(te, 1 - uu / fr)))
})

test_that("a fit carrying no orientation still errors instead of guessing", {
  d <- .a67_data()
  f <- suppressWarnings(
    sfm(y_cost ~ x1 + x2, model_name = "NHN", data = d, inefdec = FALSE)
  )
  ## A fit saved before the field existed. A43's guard has to keep holding for
  ## those, or an old cost fit reads as a production one with no warning.
  f$inefdec <- NULL
  f$call$inefdec <- as.name("no_such_variable_anywhere")
  expect_true(is.na(.is_cost_fit(f)))
  expect_error(.sfa_inefdec(f), "cannot tell whether")
  expect_error(efficiency(f, logDepVar = FALSE, newdata = d), "cannot tell whether")
})
