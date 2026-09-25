## Gap A43. Five places decided cost-vs-production orientation by four
## different mechanisms. Two were defective in the same way A41 was, and
## .sfa_inefdec() silently assumed "production" whenever it could not resolve
## the call's `inefdec` -- so a cost fit whose inefdec had gone out of scope
## was read as a production fit, with no error and no warning.

.orient_data <- function() {
  set.seed(23)
  n <- 300
  x1 <- runif(n, 1, 3)
  x2 <- runif(n, 1, 3)
  v <- rnorm(n, 0, 0.3)
  u <- abs(rnorm(n, 0, 0.6))
  data.frame(x1 = x1, x2 = x2,
    y_prod = 4 + 1.5 * x1 + 0.8 * x2 + v - u,
    y_cost = 4 + 1.5 * x1 + 0.8 * x2 + v + u)
}

.prod_fit <- function(d) suppressWarnings(sfm(y_prod ~ x1 + x2, model_name = "NHN", data = d))
.cost_fit <- function(d) {
  suppressWarnings(sfm(y_cost ~ x1 + x2, model_name = "NHN", data = d, inefdec = FALSE))
}

test_that(".sfa_inefdec() agrees with .is_cost_fit() on literal inefdec", {
  d <- .orient_data()
  fp <- .prod_fit(d)
  fc <- .cost_fit(d)
  expect_true(.sfa_inefdec(fp))
  expect_false(.sfa_inefdec(fc))
  expect_false(.is_cost_fit(fp))
  expect_true(.is_cost_fit(fc))
  ## The two helpers are inverses; that is the whole point of routing one
  ## through the other.
  expect_equal(.sfa_inefdec(fp), !.is_cost_fit(fp))
  expect_equal(.sfa_inefdec(fc), !.is_cost_fit(fc))
})

test_that("inefdec passed as a variable resolves, for a fit made inside a function", {
  d <- .orient_data()
  make <- function(dat) {
    flag <- FALSE
    suppressWarnings(sfm(y_cost ~ x1 + x2, model_name = "NHN", data = dat, inefdec = flag))
  }
  fv <- make(d)
  expect_false(.sfa_inefdec(fv))
  expect_true(.is_cost_fit(fv))
})

test_that("an unresolvable inefdec errors instead of defaulting to production", {
  d <- .orient_data()
  fx <- .cost_fit(d)
  fx$call$inefdec <- as.name("no_such_variable_anywhere")
  ## Before A43 this returned TRUE -- a cost fit silently read as production.
  expect_error(.sfa_inefdec(fx), "cannot tell whether")
  expect_error(.robust_residuals(fx), "cannot tell whether")
})

test_that(".sfa_inefdec() still honours a stored character field", {
  d <- .orient_data()
  f <- .prod_fit(d)
  f$inefdec <- "cost"
  expect_false(.sfa_inefdec(f))
  f$inefdec <- "production"
  expect_true(.sfa_inefdec(f))
})

test_that(".robust_residuals() carries the right sign for each orientation", {
  d <- .orient_data()
  fp <- .prod_fit(d)
  fc <- .cost_fit(d)
  xb <- function(f, yv) {
    b <- f$coefficients
    mm <- stats::model.matrix(~ x1 + x2, d)
    nm <- intersect(colnames(mm), names(b))
    as.numeric(d[[yv]] - mm[, nm, drop = FALSE] %*% b[nm])
  }
  expect_equal(.robust_residuals(fp), xb(fp, "y_prod"), tolerance = 1e-8)
  ## A cost frontier flips the sign: inefficiency raises y rather than
  ## lowering it, so the composed residual is -(y - Xb).
  expect_equal(.robust_residuals(fc), -xb(fc, "y_cost"), tolerance = 1e-8)
  ## The two must not be equal, or the sign branch is not being exercised.
  expect_false(isTRUE(all.equal(.robust_residuals(fc), xb(fc, "y_cost"))))
})
