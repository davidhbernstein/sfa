## psfm(derivatives = "autograd", optimizer = "newton"): exact derivatives of
## the torch-evaluated simulated likelihood, and the multistart Newton
## optimizer that uses them (gap N3). The tests pin arithmetic and invariants,
## not where an optimizer lands.

skip_if_no_torch <- function() {
  skip_on_cran()
  skip_if_not_installed("torch")
  skip_if_not(isTRUE(torch::torch_is_installed()), "libtorch not installed")
}

ad_data <- function(seed = 3) {
  data_gen_p(t = 4, N = 30, rand = seed, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
             sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
}

test_that("derivatives and optimizer are validated before any work", {
  d <- ad_data()
  expect_error(
    psfm(y_gtre ~ x1 + x2, "GTRE", d, individual = "name", estimator = "sml",
         derivatives = "autograd"),
    "needs sim_backend = \"torch\""
  )
  expect_error(
    psfm(y_gtre ~ x1 + x2, "GTRE", d, individual = "name", estimator = "sml",
         optimizer = "newton"),
    "pass derivatives = \"autograd\""
  )
  skip_if_no_torch()
  expect_error(
    psfm(y_gtre ~ x1 + x2, "GTRE", d, individual = "name", estimator = "sml",
         sim_backend = "torch", sim_device = "cpu", sim_dtype = "float32",
         derivatives = "autograd"),
    "needs double precision"
  )
})

## A cheap fit only to obtain the closure and its draws; the point at which
## the derivatives are checked is fixed by hand, not by the optimizer.
ad_closure <- function(model, inefdec = TRUE) {
  d <- ad_data()
  args <- list(formula = if (model == "GTRE") y_gtre ~ x1 + x2 else y_tre ~ x1 + x2,
               model_name = model, data = d, individual = "name",
               sim_backend = "torch", sim_device = "cpu", keep_objective = TRUE,
               optHessian = FALSE, maxit.bobyqa = 5, inefdec = inefdec)
  if (model == "GTRE") args$estimator <- "sml" else args$integration <- "simulate"
  fit <- suppressWarnings(do.call(psfm, args))
  e <- environment(fit$objective)
  ad <- sfa:::.sim_torch_autograd(e$.tprep, model, e$inefdec_n, e$n_x_vars)
  x <- if (model == "GTRE") c(2, 1, 0.3, 0.35, 0.4, 0.5, 0.5) else c(2, 1, 0.3, 0.4, 0.5, 0.5)
  list(fn = function(x) fit$objective(x), ad = ad, x = x)
}

test_that("autograd gradient and Hessian match Richardson differences, GTRE and TRE", {
  skip_if_no_torch()
  for (m in c("GTRE", "TRE")) for (dec in c(TRUE, FALSE)) {
    cl <- ad_closure(m, dec)
    g <- cl$ad$gr(cl$x)
    H <- cl$ad$he(cl$x)
    expect_equal(g, numDeriv::grad(cl$fn, cl$x), tolerance = 1e-6, info = paste(m, dec))
    expect_equal(H, numDeriv::hessian(cl$fn, cl$x), tolerance = 1e-5, info = paste(m, dec))
    expect_true(isSymmetric(H))
  }
})

test_that("the plain torch path is unchanged by the autograd refactor", {
  skip_if_no_torch()
  cl <- ad_closure("GTRE")
  e <- environment(environment(cl$fn)$fit$objective)
  r_vec <- e$fn_1_vec(cl$x, per_obs = TRUE)
  t_vec <- e$fn_1_torch(cl$x, per_obs = TRUE)
  expect_equal(t_vec, unname(r_vec), tolerance = 1e-12)
})

test_that("derivatives = 'autograd' reports the autograd Hessian at its own point", {
  skip_if_no_torch()
  d <- ad_data()
  fit <- suppressWarnings(psfm(y_gtre ~ x1 + x2, "GTRE", d, individual = "name",
    estimator = "sml", sim_backend = "torch", sim_device = "cpu",
    keep_objective = TRUE, derivatives = "autograd"))
  expect_identical(fit$derivatives, "autograd")
  expect_identical(fit$optimizer, "pipeline")
  e <- environment(fit$objective)
  ad <- sfa:::.sim_torch_autograd(e$.tprep, "GTRE", e$inefdec_n, e$n_x_vars)
  p <- unname(fit$opt$par)
  if (all(is.finite(fit$opt$hessian))) {
    expect_equal(unname(fit$opt$hessian), ad$he(p), tolerance = 1e-10)
  }
})

test_that("optimizer = 'newton' keeps the best of its starts and reports them", {
  skip_if_no_torch()
  d <- ad_data()
  for (m in c("GTRE", "TRE")) {
    args <- list(formula = if (m == "GTRE") y_gtre ~ x1 + x2 else y_tre ~ x1 + x2,
                 model_name = m, data = d, individual = "name",
                 sim_backend = "torch", sim_device = "cpu", keep_objective = TRUE,
                 derivatives = "autograd", optimizer = "newton")
    if (m == "GTRE") args$estimator <- "sml" else args$integration <- "simulate"
    fit <- suppressWarnings(do.call(psfm, args))
    st <- fit$newton$starts
    expect_identical(fit$optimizer, "newton")
    expect_equal(nrow(st), if (m == "GTRE") 6L else 3L, info = m)
    expect_equal(fit$opt$value, min(st$objective, na.rm = TRUE), info = m)
    expect_equal(st$objective[fit$newton$chosen], fit$opt$value, info = m)
    ## The reported value is the objective at the reported point.
    expect_equal(fit$objective(unname(fit$opt$par)), fit$opt$value, tolerance = 1e-10, info = m)
    ## Scales respect their floor.
    k0 <- if (m == "GTRE") 4L else 3L
    expect_true(all(fit$opt$par[seq_len(k0)] >= 1e-7), info = m)
  }
})

test_that("the Newton start grid spreads the persistent scale", {
  sv <- c(2, 1, 0.2, 0.3, 0.5, 0.5, 0.5)
  s <- sfa:::.sfa_newton_starts_panel(sv, "GTRE")
  expect_length(s, 6L)
  expect_identical(s[[1]], sv)
  ## Every GTRE alternative keeps the frontier and noise block.
  for (x in s) expect_identical(x[-c(3, 4)], sv[-c(3, 4)])
  ## A collapsed start still gets non-degenerate alternatives.
  s0 <- sfa:::.sfa_newton_starts_panel(c(2, 1, 0, 0, 0.5), "GTRE")
  expect_true(all(vapply(s0[-1], function(x) x[3] > 0 && x[4] > 0, logical(1))))
  expect_length(sfa:::.sfa_newton_starts_panel(sv[-4], "TRE"), 3L)
})
