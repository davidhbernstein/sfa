## A47. .sfa_data() used a fixed parent.frame(3) depth to re-evaluate the
## call's `data` argument. It missed the frame of a fit made inside a function,
## and -- the released, silent half of the defect -- it resolved an unrelated
## object of the same name and returned it unchecked, so fitted(), residuals(),
## predict() and efficiency(logDepVar = FALSE) all answered from the wrong data
## with no error and no warning.

.a47_data <- function(seed, cons = 5) {
  data_gen_cs(
    N = 80, rand = seed, sig_u = 0.5, sig_v = 0.2, cons = cons,
    beta1 = 1.5, beta2 = 2.0, a = 5, mu = 0.1
  )
}

## Fit on data local to this function, gone by the time it returns.
.a47_local_fit <- function() {
  d <- .a47_data(42)
  sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
}

test_that("a fit made inside a function still resolves its own data", {
  fit <- .a47_local_fit()
  truth <- .a47_data(42)
  target <- as.numeric(stats::predict(fit, newdata = truth))

  ## The frame the fit was made in is gone, but environment(object$formula)
  ## still points at it, which is what makes this recoverable at all.
  expect_equal(as.numeric(stats::fitted(fit)), target, tolerance = 1e-8)
  expect_equal(as.numeric(stats::predict(fit)), target, tolerance = 1e-8)
  expect_length(as.numeric(stats::residuals(fit)), nrow(truth))
  expect_length(as.numeric(efficiency(fit, logDepVar = FALSE)), nrow(truth))
})

test_that("a same-named decoy never displaces the fit's own data", {
  fit <- .a47_local_fit()
  truth <- .a47_data(42)
  target <- as.numeric(stats::predict(fit, newdata = truth))

  ## A different frame bound to the same name the call recorded, in the
  ## environment the old fixed-depth lookup reached. Different SEED, so x1 and
  ## x2 differ too -- a decoy that varies only the response leaves
  ## fitted() = X*beta unchanged and hides the defect entirely.
  decoy <- .a47_data(999)
  expect_false(isTRUE(all.equal(decoy$x1, truth$x1)))

  had_d <- exists("d", envir = globalenv(), inherits = FALSE)
  old_d <- if (had_d) get("d", envir = globalenv()) else NULL
  assign("d", decoy, envir = globalenv())
  on.exit(
    {
      if (had_d) assign("d", old_d, envir = globalenv()) else suppressWarnings(rm("d", envir = globalenv()))
    },
    add = TRUE
  )

  decoy_target <- as.numeric(stats::predict(fit, newdata = decoy))
  expect_false(isTRUE(all.equal(target, decoy_target, tolerance = 1e-6)))

  ## Each method must give this fit's own answer, or refuse outright -- never
  ## the decoy's answer dressed up as a success.
  fv <- tryCatch(as.numeric(stats::fitted(fit)), error = function(e) NULL)
  if (!is.null(fv)) {
    expect_equal(fv, target, tolerance = 1e-8)
    expect_false(isTRUE(all.equal(fv, decoy_target, tolerance = 1e-6)))
  }

  pv <- tryCatch(as.numeric(stats::predict(fit)), error = function(e) NULL)
  if (!is.null(pv)) {
    expect_equal(pv, target, tolerance = 1e-8)
  }

  rv <- tryCatch(as.numeric(stats::residuals(fit)), error = function(e) NULL)
  if (!is.null(rv)) {
    expect_equal(rv, as.numeric(truth$y_pcs) - target, tolerance = 1e-8)
  }

  ev <- tryCatch(as.numeric(efficiency(fit, logDepVar = FALSE)), error = function(e) NULL)
  if (!is.null(ev)) {
    expect_equal(ev, as.numeric(efficiency(fit, newdata = truth, logDepVar = FALSE)),
      tolerance = 1e-8
    )
  }
})

test_that("an unrecoverable fit still errors clearly rather than succeeding", {
  ## Nothing named `d` anywhere, and the formula's environment emptied, so
  ## there is genuinely nothing to find. The fix must not turn this clear
  ## message into a silent success against some other frame.
  fit <- .a47_local_fit()
  environment(fit$formula) <- new.env(parent = emptyenv())

  had_d <- exists("d", envir = globalenv(), inherits = FALSE)
  old_d <- if (had_d) get("d", envir = globalenv()) else NULL
  if (had_d) rm("d", envir = globalenv())
  on.exit(if (had_d) assign("d", old_d, envir = globalenv()), add = TRUE)

  expect_error(stats::fitted(fit), "Cannot recover the data")
  expect_error(stats::predict(fit), "Cannot recover the data")
  expect_error(stats::residuals(fit), "Cannot recover the data")
  expect_error(efficiency(fit, logDepVar = FALSE), "could not rebuild the fitted frontier")
})

test_that("nobs() does not take its count from a same-named decoy", {
  fit <- .a47_local_fit()

  ## Force the call-re-evaluation branch: strip every field nobs.sfareg()
  ## would otherwise answer from, so the environment search is what is under
  ## test rather than a stored count.
  stripped <- fit
  stripped$nobs <- NULL
  stripped$data <- NULL
  stripped$exp_u_hat <- NULL
  stripped$u_hat <- NULL
  stripped$residuals <- NULL
  stripped$med_u_hat <- NULL
  stripped$u_posterior <- NULL

  decoy <- .a47_data(999)[1:37, , drop = FALSE]
  had_d <- exists("d", envir = globalenv(), inherits = FALSE)
  old_d <- if (had_d) get("d", envir = globalenv()) else NULL
  assign("d", decoy, envir = globalenv())
  on.exit(
    {
      if (had_d) assign("d", old_d, envir = globalenv()) else suppressWarnings(rm("d", envir = globalenv()))
    },
    add = TRUE
  )

  n <- stats::nobs(stripped)
  expect_false(identical(as.integer(n), 37L))
  expect_true(is.na(n) || identical(as.integer(n), 80L))
})

## A44. The fit drops rows on complete.cases() across every pipe segment at
## once, response included (data.processing.R:219). The rebuild in .sfa_xb()
## saw only the frontier, and delete.response() hid y from it, so rows the fit
## never used survived into fitted()/residuals()/predict() -- one value too
## many, silently misaligned against exp_u_hat. residuals() went further and
## recycled the response against a shorter design, warning but returning
## nonsense.

test_that("rebuilt vectors use the rows the fit used, not the rows supplied", {
  base <- data_gen_cs(
    N = 120, rand = 11, sig_u = 0.5, sig_v = 0.2, cons = 5,
    beta1 = 1.5, beta2 = 2.0, a = 5, mu = 0.1
  )[, c("y_pcs_z", "x1", "x2", "z")]

  ## One NA at a time, in each role a variable can play: a variance
  ## determinant outside rhs = 1, the response that delete.response() hides,
  ## and an ordinary frontier regressor.
  for (col in c("z", "y_pcs_z", "x1")) {
    d <- base
    d[[col]][7] <- NA

    fit <- suppressWarnings(
      sfm(y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z", data = d)
    )
    n_used <- length(as.numeric(fit$exp_u_hat))
    expect_equal(n_used, 119L, info = col)

    expect_length(as.numeric(stats::fitted(fit)), n_used)
    expect_length(as.numeric(stats::predict(fit)), n_used)
    expect_length(as.numeric(stats::residuals(fit)), n_used)
    ## This design puts a few fitted frontiers at or below zero, which
    ## efficiency() warns about on the level scale. Expected here and not what
    ## this test is about.
    expect_length(
      suppressWarnings(as.numeric(efficiency(fit, logDepVar = FALSE))), n_used
    )

    ## residuals() recycled rather than erroring, so a length check alone is
    ## not enough -- the values must be the real ones.
    expect_silent(rv <- as.numeric(stats::residuals(fit)))
    expect_true(all(is.finite(rv)), info = col)
  }
})
