## efficiency_ci()'s te_hat fallback, for a fit that no longer carries
## exp_u_hat.
##
## The reference is the Battese-Coelli definition integrated numerically:
## E[exp(-u) | e] for u | e truncated normal on (0, Inf). That is an
## independent implementation of the quantity rather than the fit's own stored
## field, so these tests pin the estimator instead of pinning agreement with
## whatever the fit happened to save. Closed form and quadrature agree to
## 2e-12 over mu_star in [-2, 5].
##
## Before A68 the fallback computed exp(-JLMS) = exp(-E[u | e]), which Jensen
## separates from E[exp(-u) | e] -- a measured 8.0e-3 on the NHN fixture here.

## E[exp(-u)] for u ~ N(mu, s^2) truncated below at zero. Plain-scale
## quadrature, usable only where integrate() can actually resolve the density:
## NHN_Z collapses sigma_star to 5.3e-07 on this fixture (that is A63), and a
## spike that narrow integrates to zero over (0, Inf). The package's own closed
## form works in logs and has no such limit, which is exactly why it cannot be
## its own reference.
bc_quadrature <- function(mu, s) {
  vapply(seq_along(mu), function(i) {
    num <- stats::integrate(
      function(u) exp(-u) * stats::dnorm(u, mu[i], s[i]),
      lower = 0, upper = Inf, rel.tol = 1e-11
    )$value
    num / stats::pnorm(mu[i] / s[i])
  }, numeric(1))
}

## Observations where that quadrature is trustworthy.
a68_reliable <- function(post, n = 5L) {
  mu <- as.numeric(post$mu_star)
  s <- as.numeric(post$sigma_star)
  ok <- which(is.finite(mu / s) & mu / s > -5 & s > 1e-6)
  utils::head(ok, n)
}

a68_fits <- function() {
  d <- cs_small(N = 150)
  ## NHN_Z stops at a saddle on this fixture -- A63, not A68. The fallback has
  ## to behave whatever the optimizer did, so the fit is kept and the warning
  ## silenced here rather than the model dropped.
  suppressWarnings(list(
    NHN = sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d),
    NHN_Z = sfm(y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z", data = d),
    NE = sfm(y_pcs ~ x1 + x2, model_name = "NE", data = d),
    NTN = sfm(y_pcs_tn ~ x1 + x2, model_name = "NTN", data = d)
  ))
}

test_that("the te_hat fallback reports Battese-Coelli, not exp(-JLMS)", {
  skip_on_cran()
  fits <- a68_fits()
  checked <- character(0)
  for (nm in names(fits)) {
    post <- fits[[nm]]$u_posterior
    ## A fit saved before $exp_u_hat existed, or stripped of it.
    stripped <- fits[[nm]]
    stripped$exp_u_hat <- NULL
    ci <- efficiency_ci(stripped)

    ## Shape and range hold for every model, reference or no reference.
    expect_length(ci$te_hat, length(as.numeric(post$mu_star)))
    expect_false(is.matrix(ci$te_hat), info = nm)
    expect_true(all(ci$te_hat >= 0 & ci$te_hat <= 1), info = nm)

    i <- a68_reliable(post)
    if (length(i) < 2L) next
    expect_equal(ci$te_hat[i],
                 bc_quadrature(as.numeric(post$mu_star)[i],
                               as.numeric(post$sigma_star)[i]),
                 tolerance = 1e-8, info = nm)
    checked <- c(checked, nm)
  }
  ## Guard against the reference quietly going unused: say which models the
  ## definition was actually checked against.
  expect_setequal(checked, c("NHN", "NE", "NTN"))
})

test_that("the fallback and the stored field cannot drift apart", {
  skip_on_cran()
  ## The no-drift guarantee the function's own comment claims: what
  ## efficiency_ci() reports must not depend on whether exp_u_hat survived.
  ## Exact rather than approximate, because the fit built that field from this
  ## same posterior -- and this is the one block NHN_Z is covered by.
  fits <- a68_fits()
  for (nm in names(fits)) {
    stripped <- fits[[nm]]
    stripped$exp_u_hat <- NULL
    expect_equal(efficiency_ci(stripped)$te_hat,
                 efficiency_ci(fits[[nm]])$te_hat,
                 tolerance = 1e-12, info = nm)
  }
})

test_that("the fixture separates the two estimators, so the check discriminates", {
  skip_on_cran()
  ## Without this the blocks above would pass on the old code wherever the fit
  ## landed somewhere the two estimators coincide -- which they do as
  ## sigma_star -> 0, NHN_Z's collapse here being exactly such a place. Pin
  ## that NHN is not such a place rather than assuming it.
  fit <- a68_fits()$NHN
  post <- fit$u_posterior
  jlms <- exp(-sfa:::.jlms_u(post$mu_star, post$sigma_star))
  bc <- as.numeric(fit$exp_u_hat)
  gap <- max(abs(bc - jlms))
  skip_if_not(gap > 1e-3,
              sprintf("fit landed where the two estimators agree (gap %.2g)", gap))

  stripped <- fit
  stripped$exp_u_hat <- NULL
  te <- efficiency_ci(stripped)$te_hat
  ## Jensen: exp(-E[u]) <= E[exp(-u)], so the old fallback was biased DOWN.
  expect_true(all(jlms <= bc + 1e-12))
  expect_gt(max(abs(te - jlms)), 1e-3)
  expect_equal(te, bc, tolerance = 1e-12)
})
