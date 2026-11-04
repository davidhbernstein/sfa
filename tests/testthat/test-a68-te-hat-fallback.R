## efficiency_ci()'s te_hat fallback, for a fit that no longer carries
## exp_u_hat.
##
## Before A68 the fallback computed exp(-JLMS) = exp(-E[u | e]) where the
## stored field is the Battese-Coelli score E[exp(-u) | e]. Jensen separates
## them, so the column reported a different estimator under the same name.
##
## The estimator is pinned on a posterior written down here rather than on a
## fitted one. That is deliberate: it makes the reference exact and the test
## platform-independent. An earlier version of this file drove the same check
## from fits and selected the observations where plain-scale quadrature could
## be trusted. It passed locally and failed on four of five CI platforms,
## because which observations those are depends on where the optimizer landed
## -- NHN_Z's sigma_star collapses to 5.3e-07 here (that is A63) and to
## something larger on CI. The fits below are still used, but only for
## invariants that hold wherever they land.

## E[exp(-u)] for u ~ N(mu, s^2) truncated below at zero, by quadrature. NA
## when integrate() cannot resolve the density -- for a narrow spike far from
## zero it silently returns 0, which is why the convergence message is checked
## rather than the value trusted.
bc_quadrature <- function(mu, s) {
  vapply(seq_along(mu), function(i) {
    r <- try(stats::integrate(
      function(u) exp(-u) * stats::dnorm(u, mu[i], s[i]),
      lower = 0, upper = Inf, rel.tol = 1e-11
    ), silent = TRUE)
    if (inherits(r, "try-error") || !identical(r$message, "OK")) {
      return(NA_real_)
    }
    r$value / stats::pnorm(mu[i] / s[i])
  }, numeric(1))
}

test_that("the te_hat fallback is Battese-Coelli, not exp(-JLMS)", {
  ## A posterior chosen so the quadrature is reliable and the two estimators
  ## are far apart. efficiency_ci() needs nothing else from the object, so no
  ## optimizer is involved and the expected values are the same everywhere.
  mu <- c(-1.0, -0.4, 0.0, 0.3, 0.8, 1.5, 2.0)
  s <- c(0.5, 0.7, 0.4, 0.9, 0.6, 1.0, 0.8)
  obj <- structure(
    list(model_name = "NHN",
         u_posterior = list(mu_star = mu, sigma_star = s)),
    class = "sfareg"
  )
  te <- efficiency_ci(obj)$te_hat

  ## Against the definition.
  want <- bc_quadrature(mu, s)
  expect_false(anyNA(want))
  expect_equal(te, want, tolerance = 1e-10)

  ## And not the other estimator. Jensen fixes the direction: exp(-E[u]) is
  ## below E[exp(-u)], so the old fallback was biased DOWN at every point.
  jlms <- exp(-sfa:::.jlms_u(mu, s))
  expect_true(all(jlms < te))
  expect_gt(min(te - jlms), 1e-2)

  ## A probability per observation, and not a one-column matrix -- mu_star
  ## arrives as a matrix from some models and would otherwise change the
  ## column's type.
  expect_length(te, length(mu))
  expect_false(is.matrix(te))
  expect_true(all(te >= 0 & te <= 1))
})

test_that("what efficiency_ci reports does not depend on exp_u_hat surviving", {
  skip_on_cran()
  ## The no-drift guarantee, on real fits, for every model with a
  ## truncated-normal posterior. Exact, because the fit built that field from
  ## this same posterior -- so this holds wherever the optimizer landed, which
  ## is what makes it safe to assert on every platform.
  d <- cs_small(N = 150)
  ## NHN_Z stops at a saddle on this fixture (A63, not this gap). The fallback
  ## has to behave whatever the optimizer did, so the fit is kept.
  fits <- suppressWarnings(list(
    NHN = sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d),
    NHN_Z = sfm(y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z", data = d),
    NE = sfm(y_pcs ~ x1 + x2, model_name = "NE", data = d),
    NTN = sfm(y_pcs_tn ~ x1 + x2, model_name = "NTN", data = d)
  ))
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    stripped <- fit
    stripped$exp_u_hat <- NULL
    te <- efficiency_ci(stripped)$te_hat

    expect_equal(te, efficiency_ci(fit)$te_hat, tolerance = 1e-12, info = nm)
    expect_length(te, length(as.numeric(fit$u_posterior$mu_star)))
    expect_false(is.matrix(te), info = nm)
    expect_true(all(te >= 0 & te <= 1), info = nm)
    ## The stored field is genuinely present, so the comparison above is not a
    ## vacuous NULL-against-NULL.
    expect_false(is.null(fit$exp_u_hat), info = nm)
  }
})
