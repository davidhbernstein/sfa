## Gap M4. marginal_effects() covered the families whose E[u] is proportional to
## the scale, where the effect is half * delta * E[u]. The truncated normal is
## not one: mu enters separately, so the derivative is genuinely two-term. These
## tests check the ALGEBRA against numerical differentiation of the closed forms,
## which is where a hand-derived chain rule actually goes wrong.

## The closed forms, written out independently of the package's helpers.
E_tn <- function(mu, sig) {
  a <- mu / sig
  mu + sig * (dnorm(a) / pnorm(a))
}
V_tn <- function(mu, sig) {
  a <- mu / sig
  lam <- dnorm(a) / pnorm(a)
  sig^2 * (1 - lam * (lam + a))
}

test_that("dE/dmu and dE/dsigma match numerical differentiation", {
  grid <- expand.grid(mu = c(-1.5, -0.4, 0, 0.3, 1.2, 3), sig = c(0.2, 0.7, 1, 2.5))
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]
    sig <- grid$sig[i]
    a <- mu / sig
    lam <- .sfa_lambda_ratio(a)
    g <- 1 - lam * (lam + a)

    analytic_dmu <- g
    analytic_dsig <- lam + a * lam * (lam + a)

    num_dmu <- numDeriv::grad(function(m) E_tn(m, sig), mu)
    num_dsig <- numDeriv::grad(function(s) E_tn(mu, s), sig)

    expect_lt(abs(analytic_dmu - num_dmu), 1e-6 * max(1, abs(num_dmu)))
    expect_lt(abs(analytic_dsig - num_dsig), 1e-6 * max(1, abs(num_dsig)))
  }
})

test_that("g'(a), the Var factor's derivative, matches numerical differentiation", {
  for (a in c(-1.5, -0.4, 0, 0.3, 1.2, 3)) {
    lam <- .sfa_lambda_ratio(a)
    gp_analytic <- lam * (lam + a) * (2 * lam + a) - lam
    g_of <- function(aa) {
      l <- dnorm(aa) / pnorm(aa)
      1 - l * (l + aa)
    }
    gp_num <- numDeriv::grad(g_of, a)
    expect_lt(abs(gp_analytic - gp_num), 1e-6 * max(1, abs(gp_num)))
  }
})

test_that(".sfa_lambda_ratio() is the inverse Mills ratio, and survives the left tail", {
  for (a in c(-6, -3, -1, 0, 1, 4)) {
    expect_lt(abs(.sfa_lambda_ratio(a) - dnorm(a) / pnorm(a)), 1e-10)
  }
  ## Where the naive ratio underflows to 0/0 = NaN, the log form must not.
  expect_true(is.finite(.sfa_lambda_ratio(-40)))
  expect_gt(.sfa_lambda_ratio(-40), 0)
})

test_that(".sfa_me_delta_by_name() matches by name and zero-fills the absent", {
  spec <- list(delta = c(z1 = 0.5, z3 = -0.2))
  got <- .sfa_me_delta_by_name(spec, c("z1", "z2", "z3"))
  expect_equal(unname(got), c(0.5, 0, -0.2))
  expect_identical(names(got), c("z1", "z2", "z3"))
  ## A covariate in NO design contributes zero, it does not error.
  expect_equal(unname(.sfa_me_delta_by_name(NULL, c("z1"))), 0)
  expect_equal(unname(.sfa_me_delta_by_name(list(), c("z1"))), 0)
})

test_that("a truncated-normal het fit reports both terms, matching a hand chain rule", {
  skip_on_cran()
  set.seed(404)
  n <- 250
  d <- data_gen_cs(N = n, rand = 404, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  d$z <- runif(n, -1, 1)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NTN", data = d, uhet = ~z, muhet = ~z)
  me <- marginal_effects(f)
  expect_s3_class(me, "data.frame")
  expect_identical(attr(me, "family"), "truncnormal")
  expect_true(any(grepl("^dE_u[.]dz$", names(me))))

  ## Rebuild the effect from the fit's own pieces, independently of the code
  ## under test: sigma_u and mu from the stored designs, then numeric gradients.
  zs <- f$z_spec
  ms <- f$mu_spec
  half <- if (identical(zs$link, "sd")) 1 else 0.5
  sig_of <- function(zv, row) {
    Zr <- zs$Z[row, , drop = FALSE]; Zr[, "z"] <- zv
    eta <- as.numeric(Zr %*% zs$delta)
    if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))
  }
  mu_of <- function(zv, row) {
    Zr <- ms$Z[row, , drop = FALSE]; Zr[, "z"] <- zv
    as.numeric(Zr %*% ms$delta)
  }
  for (row in c(1L, 50L, 200L)) {
    num <- numDeriv::grad(function(zv) E_tn(mu_of(zv, row), sig_of(zv, row)),
                          zs$Z[row, "z"])
    expect_lt(abs(me[[grep("^dE_u[.]dz$", names(me))]][row] - num),
              1e-5 * max(1, abs(num)))
  }
})

test_that("a covariate in the mu design only is reported, not dropped", {
  skip_on_cran()
  ## The covariate list used to come from the sigma_u design alone, so with
  ## muhet = ~ z2 and uhet = ~ z1 the z2 effect was silently missing, and with
  ## muhet alone the call stopped claiming every column was constant.
  set.seed(405)
  n <- 250
  d <- data_gen_cs(N = n, rand = 405, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  d$z1 <- runif(n, -1, 1)
  d$z2 <- runif(n, -1, 1)

  f <- sfm(y_pcs ~ x1 + x2, model_name = "NTN", data = d,
           uhet = ~z1, muhet = ~z2)
  me <- marginal_effects(f)
  expect_true(all(c("dE_u.dz1", "dE_u.dz2", "dVar_u.dz1", "dVar_u.dz2") %in%
                    names(me)))

  ## Independent check: numerical derivatives of the closed-form moments,
  ## moving one covariate in its own design.
  zs <- f$z_spec
  ms <- f$mu_spec
  sig_of <- function(z1, row) {
    Zr <- zs$Z[row, , drop = FALSE]; Zr[, "z1"] <- z1
    eta <- as.numeric(Zr %*% zs$delta)
    if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))
  }
  mu_of <- function(z2, row) {
    Zr <- ms$Z[row, , drop = FALSE]; Zr[, "z2"] <- z2
    as.numeric(Zr %*% ms$delta)
  }
  for (row in c(1L, 50L, 200L)) {
    s0 <- sig_of(zs$Z[row, "z1"], row)
    m0 <- mu_of(ms$Z[row, "z2"], row)
    nE1 <- numDeriv::grad(function(v) E_tn(m0, sig_of(v, row)), zs$Z[row, "z1"])
    nE2 <- numDeriv::grad(function(v) E_tn(mu_of(v, row), s0), ms$Z[row, "z2"])
    nV2 <- numDeriv::grad(function(v) V_tn(mu_of(v, row), s0), ms$Z[row, "z2"])
    expect_lt(abs(me$dE_u.dz1[row] - nE1), 1e-5 * max(1, abs(nE1)))
    expect_lt(abs(me$dE_u.dz2[row] - nE2), 1e-5 * max(1, abs(nE2)))
    expect_lt(abs(me$dVar_u.dz2[row] - nV2), 1e-5 * max(1, abs(nV2)))
  }

  ## mu-only heterogeneity: sigma_u is constant, the effect runs through mu.
  g <- sfm(y_pcs ~ x1 + x2, model_name = "NTN", data = d, muhet = ~z2)
  me2 <- marginal_effects(g)
  expect_true("dE_u.dz2" %in% names(me2))
  expect_true(all(is.finite(me2$dE_u.dz2)))
})
