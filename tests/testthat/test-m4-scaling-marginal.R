## Gap M4, item 1. Under sfm(scaling = ~ z) one factor h = exp(z'delta_s)
## multiplies BOTH the scale and the pre-truncation mean, so a = mu/sigma_u is
## constant in z, E[u] = h E[u*], and the marginal effect is exactly
## delta_s_k * E[u] -- with delta_s_k itself the semi-elasticity. The fit stored
## no scaling block, which is the only reason this could not be reported.

make_scaling_fit <- function(n = 300, seed = 9) {
  d <- data_gen_cs(N = n, rand = seed, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  set.seed(seed)
  d$z <- runif(n, -1, 1)
  ## Scaling is implemented for NTN only: for a half-normal u the constraint
  ## adds nothing, since h*|N(0,s^2)| IS |N(0,(h s)^2)|.
  sfm(y_pcs ~ x1 + x2, model_name = "NTN", data = d, scaling = ~z)
}

test_that("the fit stores a scaling block whose deltas ARE its scale.* rows", {
  skip_on_cran()
  f <- make_scaling_fit()
  expect_false(is.null(f$s_spec))
  sc <- grep("^scale\\.", rownames(f$out), value = TRUE)
  expect_length(sc, 1L)
  ## The offset into the parameter vector is derived, not given by n_blocks, so
  ## it is pinned against the fit's own named rows rather than assumed.
  expect_equal(unname(f$s_spec$delta), unname(f$out[sc, "par"]))
  expect_identical(names(f$s_spec$delta), "z")
})

test_that("the scaling marginal effect is exactly delta_s * E[u]", {
  skip_on_cran()
  f <- make_scaling_fit()
  me <- marginal_effects(f)
  ds <- unname(f$s_spec$delta[["z"]])
  eu <- me[[grep("^E_u$", names(me))]]
  got <- me[[grep("^dE_u[.]dz$", names(me))]]
  expect_equal(got, ds * eu, tolerance = 1e-12)
  gotv <- me[[grep("^dVar_u[.]dz$", names(me))]]
  expect_equal(gotv, 2 * ds * me[[grep("^Var_u$", names(me))]], tolerance = 1e-12)
})

test_that("it agrees with numerical differentiation of E[u] through h", {
  skip_on_cran()
  f <- make_scaling_fit()
  me <- marginal_effects(f)
  ss <- f$s_spec
  zs <- f$z_spec
  ms <- f$mu_spec
  s0 <- if (identical(zs$link, "sd")) exp(as.numeric(zs$Z[1, ] %*% zs$delta)) else
    sqrt(exp(as.numeric(zs$Z[1, ] %*% zs$delta)))
  m0 <- as.numeric(ms$Z[1, ] %*% ms$delta)
  ## E[u](z) built from the definition, with h applied to both scale and mean.
  Eu_of <- function(zv) {
    h <- exp(zv * unname(ss$delta[["z"]]))
    mu <- m0 * h; sig <- s0 * h; a <- mu / sig
    mu + sig * (dnorm(a) / pnorm(a))
  }
  col <- grep("^dE_u[.]dz$", names(me))
  for (row in c(1L, 77L, 200L)) {
    num <- numDeriv::grad(Eu_of, ss$Z[row, "z"])
    expect_lt(abs(me[[col]][row] - num), 1e-6 * max(1, abs(num)))
  }
})

test_that("delta_s is the semi-elasticity, i.e. dE[u]/dz / E[u] is constant", {
  skip_on_cran()
  f <- make_scaling_fit()
  me <- marginal_effects(f)
  r <- me[[grep("^dE_u[.]dz$", names(me))]] / me[[grep("^E_u$", names(me))]]
  expect_lt(diff(range(r)), 1e-12)
  expect_equal(unname(r[1]), unname(f$s_spec$delta[["z"]]), tolerance = 1e-10)
})

test_that("a non-scaling fit is untouched by the new branch", {
  skip_on_cran()
  d <- data_gen_cs(N = 250, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  set.seed(5); d$z <- runif(250, -1, 1)
  f <- sfm(y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z", data = d)
  expect_null(f$s_spec)
  me <- marginal_effects(f)
  half <- if (identical(f$z_spec$link, "sd")) 1 else 0.5
  expect_equal(me[[grep("^dE_u[.]dz$", names(me))]],
               half * unname(f$z_spec$delta[["z"]]) * me[[grep("^E_u$", names(me))]],
               tolerance = 1e-12)
})
