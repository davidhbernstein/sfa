## psfm(sim_backend = "torch"): the optional torch evaluator of the GTRE (sml)
## and TRE simulated likelihoods.

test_that("the torch arguments are validated before any work", {
  d <- data_gen_p(t = 4, N = 20, rand = 3, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  ## Not one of the two stacked simulated likelihoods.
  expect_error(
    suppressWarnings(psfm(y_gtre ~ x1 + x2, "GTRE_Z", pd, individual = "name",
                          sim_backend = "torch")),
    "implemented for model_name"
  )
  expect_error(
    suppressWarnings(psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name",
                          estimator = "fiml", sim_backend = "torch")),
    "implemented for model_name"
  )
  ## Device and precision only mean something with torch.
  expect_warning(
    sfa:::.sim_torch_setup("R", "cpu", NULL, "GTRE"),
    "apply only with"
  )
  expect_null(suppressWarnings(sfa:::.sim_torch_setup("R", NULL, NULL, "GTRE")))
})

skip_if_no_torch <- function() {
  skip_if_not_installed("torch")
  skip_if_not(isTRUE(torch::torch_is_installed()), "libtorch not installed")
}

test_that("the torch log-density matches .gtre_sim_logdens() exactly in double", {
  skip_if_no_torch()
  set.seed(1)
  n <- 40; R <- 7; N <- 8
  E <- matrix(rnorm(n * R, sd = 2), n, R)
  gid <- rep(seq_len(N), each = n / N)
  st <- sfa:::.sim_torch_setup("torch", "cpu", "float64", "GTRE")
  prep <- list(setup = st, N = N, R = R,
               gid = torch::torch_tensor(gid, dtype = torch::torch_long()))
  r_ll <- sfa:::.gtre_sim_logdens(E, 1.7, 0.9, gid, N)
  t_ll <- as.numeric(sfa:::.gtre_sim_logdens_torch(
    torch::torch_tensor(E, dtype = torch::torch_float64()), 1.7, 0.9, prep))
  ## A leading R scalar in a torch expression is converted to float32; this
  ## is the test that catches it (it showed up as ~1e-8).
  expect_equal(t_ll, unname(r_ll), tolerance = 1e-13)
})

test_that("single-precision log Phi is accurate through the clip range", {
  skip_if_no_torch()
  z <- seq(-37, 8, by = 0.05)
  t32 <- as.numeric(sfa:::.torch_log_pnorm(
    torch::torch_tensor(z, dtype = torch::torch_float32()), single = TRUE))
  ref <- pnorm(z, log.p = TRUE)
  ## Relative to |log Phi| where it is large, absolute near zero.
  expect_true(all(is.finite(t32)))
  expect_lt(max(abs(t32 - ref) / pmax(1, abs(ref))), 1e-5)
})

test_that("per-firm log-likelihoods agree across backends, both orientations", {
  skip_if_no_torch()
  skip_on_cran()
  d <- data_gen_p(t = 6, N = 40, rand = 7, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  fit_obj <- function(m, prod, ...) {
    pd <- plm::pdata.frame(d, index = c("name", "year"))
    if (!prod) { pd$y_gtre <- -pd$y_gtre; pd$y_tre <- -pd$y_tre }
    f <- if (m == "GTRE") y_gtre ~ x1 + x2 else y_tre ~ x1 + x2
    a <- list(f, m, pd, individual = "name", halton_num = 30, rand.gtre = 7L,
              keep_objective = TRUE, inefdec = prod, maxit.bobyqa = 40,
              maxit.optim = 1, optHessian = FALSE, ...)
    if (m == "GTRE") a$estimator <- "sml"
    suppressWarnings(do.call(psfm, a))
  }
  for (m in c("GTRE", "TRE")) for (prod in c(TRUE, FALSE)) {
    a <- fit_obj(m, prod)
    x <- a$opt$par
    ref <- a$objective(x, per_obs = TRUE)
    t64 <- fit_obj(m, prod, sim_backend = "torch", sim_device = "cpu")
    t32 <- fit_obj(m, prod, sim_backend = "torch", sim_device = "cpu",
                   sim_dtype = "float32")
    expect_equal(t64$objective(x, per_obs = TRUE), ref, tolerance = 1e-12,
                 info = paste(m, prod))
    expect_lt(max(abs(t32$objective(x, per_obs = TRUE) - ref)), 1e-4)
    expect_identical(t64$sim_backend,
                     list(backend = "torch", device = "cpu", dtype = "float64"))
    expect_identical(a$sim_backend$backend, "R")
  }
})

test_that("a full double-precision torch fit reproduces the default fit", {
  skip_if_no_torch()
  skip_on_cran()
  d <- data_gen_p(t = 5, N = 50, rand = 11, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  ## Same RNG state for both, so psoptim's swarm is the same in each.
  set.seed(42)
  a <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pd, individual = "name",
                             halton_num = 30, rand.gtre = 7L))
  set.seed(42)
  b <- suppressWarnings(psfm(y_tre ~ x1 + x2, "TRE", pd, individual = "name",
                             halton_num = 30, rand.gtre = 7L,
                             sim_backend = "torch", sim_device = "cpu"))
  ## Same objective to rounding error: the invariant is the log-likelihood at
  ## each one's optimum, not where the optimizer lands.
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-6)
})
