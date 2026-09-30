## Gap A52/A52a: optHessian = FALSE, PSopt = FALSE.
##
## Without the final optim() stage, a fit stored whichever earlier stage ran as
## `$opt`, raw. Two of the four stages do not use optim()'s field names --
## nlminb() returns `objective`, bobyqa() returns `fval`/`ierr` -- and every
## consumer reads `$value`. Measured on main before the fix:
##
##   sfm()   NHN NE NTN NU   hard error: nlminb ran, its result was discarded,
##                           `opt` was NULL and `out[1, ] <- opt$par` failed
##   sfm()   the other 11    logLik() silently NA (bobyqa object stored)
##   lcsfm() LCM LCM_Z       hard error: `-opt$value` on NULL
##   lcsfm() LCM_CN          hard error: `opt` never set at all

cs_a52 <- function() cs_small(N = 200, rand = 1)

test_that(".as_optim() gives optim()'s names to nlminb and bobyqa results", {
  fn <- function(p) sum((p - c(1, 2))^2) + 3
  nl <- stats::nlminb(c(0, 0), fn)
  bq <- minqa::bobyqa(c(0, 0), fn, lower = c(-5, -5), upper = c(5, 5),
    control = list(rhobeg = 0.5, rhoend = 1e-8))

  a <- .as_optim(nl)
  expect_identical(a$value, nl$objective)
  expect_identical(a$par, nl$par)
  expect_identical(a$convergence, as.integer(nl$convergence))

  b <- .as_optim(bq)
  expect_identical(b$value, bq$fval)
  expect_identical(b$par, bq$par)
  expect_identical(b$convergence, as.integer(bq$ierr))

  ## Both are the objective at the point they report -- checked against the
  ## objective itself, not against another field of the same object.
  expect_equal(a$value, fn(a$par))
  expect_equal(b$value, fn(b$par))

  ## optim() and psoptim() already conform and must pass through untouched.
  op <- stats::optim(c(0, 0), fn, method = "BFGS", hessian = TRUE)
  expect_identical(.as_optim(op), op)
  expect_null(.as_optim(NULL))
})

test_that("sfm()'s nlminb models fit, and report the likelihood at their estimate", {
  skip_on_cran()
  d <- cs_a52()
  for (m in c("NHN", "NE", "NTN", "NU")) {
    y <- switch(m, NHN = "y_pcs", NE = "y_pcs_e", NTN = "y_pcs_tn", NU = "y_pcs_u")
    f <- sfm(stats::as.formula(paste(y, "~ x1 + x2")), model_name = m, data = d,
      optHessian = FALSE)
    expect_true(is.finite(as.numeric(logLik(f))), info = m)
    expect_true(all(is.finite(f$out[, "par"])), info = m)
  }

  ## NHN against the closed-form density, at the reported estimate:
  ## f(e) = (2 / sigma) phi(e / sigma) Phi(-lambda e / sigma).
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, optHessian = FALSE)
  p <- f$out[, "par"]
  e <- d$y_pcs - cbind(1, d$x1, d$x2) %*% p[3:5]
  ll <- sum(log(2) - log(p[2]) + dnorm(e / p[2], log = TRUE) +
    pnorm(-p[1] * e / p[2], log.p = TRUE))
  expect_equal(as.numeric(logLik(f)), ll, tolerance = 1e-8)
})

test_that("sfm()'s bobyqa models no longer lose the log-likelihood", {
  skip_on_cran()
  d <- cs_a52()
  f <- sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d, optHessian = FALSE)
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_false(is.null(f$opt$convergence))

  ## And it is the likelihood at the reported estimate, through the fit's own
  ## objective: -value is the log-likelihood of the point in out[, "par"].
  f2 <- sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d, optHessian = FALSE,
    keep_objective = TRUE)
  expect_equal(as.numeric(logLik(f2)), -f2$objective(f2$out[, "par"]),
    tolerance = 1e-8)
})

test_that("lcsfm() fits all three models without the final stage", {
  skip_on_cran()
  set.seed(21)
  n <- 300
  x1 <- rnorm(n)
  cls <- rbinom(n, 1, 0.5)
  d <- data.frame(
    y = ifelse(cls == 1, 4, 1) + x1 + rnorm(n, 0, 0.4) - abs(rnorm(n, 0, 0.6)),
    x1 = x1
  )
  for (m in c("LCM", "LCM_Z", "LCM_CN")) {
    frm <- if (m == "LCM_Z") y ~ x1 | x1 else y ~ x1
    f <- lcsfm(frm, model_name = m, data = d, n_class = 2, optHessian = FALSE)
    expect_true(is.finite(as.numeric(logLik(f))), info = m)
  }
})
