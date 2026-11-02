## efficiency_ci()'s guards on the predictors it takes from the fit.
##
## Two guards, neither of which had a test before A69. The length guard was
## added for GitHub issue #2, where a column of the wrong length is recycled by
## the data frame and misaligned silently. The identity guard is A69: length
## alone accepts a vector of the right length and the wrong order, which is
## still wrong and is not caught by anything downstream.
##
## Each stored predictor is a fixed function of the posterior the fit also
## carries, so the identity is exact and checkable. Swept before relying on it:
## 72 fits (N in 80/150/300, seeds 1/42/99, NHN/NHN_Z/NE/NTN, both
## orientations) put the worst discrepancy at 7.77e-16 for exp_u_hat and
## exactly 0 for u_hat, against a guard tolerance of 1e-8. So the guard has
## three orders of magnitude of headroom on real fits and still catches
## corruption at the 1e-2 scale a reordering produces.

a69_fit <- function(N = 150) {
  sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = N))
}

## A posterior written down rather than fitted, so the guard's behaviour is
## pinned identically on every platform. The fitted cases below add realism;
## this one adds certainty.
a69_handmade <- function() {
  mu <- c(-1.0, -0.4, 0.0, 0.3, 0.8, 1.5, 2.0)
  s <- c(0.5, 0.7, 0.4, 0.9, 0.6, 1.0, 0.8)
  obj <- structure(
    list(model_name = "NHN",
         u_posterior = list(mu_star = mu, sigma_star = s)),
    class = "sfareg"
  )
  obj$exp_u_hat <- as.numeric(sfa:::.te_battese_coelli(mu, s))
  obj$u_hat <- as.numeric(sfa:::.jlms_u(mu, s))
  obj
}

test_that("the guard is exact on a posterior written down, not fitted", {
  obj <- a69_handmade()
  ## Correct by construction, so it must pass.
  expect_s3_class(efficiency_ci(obj), "data.frame")

  ## Reversal: same values, same mean, same range, pairing destroyed.
  rev_obj <- obj
  rev_obj$exp_u_hat <- rev(obj$exp_u_hat)
  expect_error(efficiency_ci(rev_obj), "exp_u_hat")
  expect_error(efficiency_ci(rev_obj), "do not belong to these observations")

  ## One observation moved by 1e-4, against a guard tolerance of 1e-8. Both
  ## numbers are fixed here, so this sensitivity does not depend on a fit.
  one <- obj
  one$exp_u_hat[3] <- one$exp_u_hat[3] + 1e-4
  expect_error(efficiency_ci(one), "do not belong to these observations")

  ## And the same for u_hat.
  u_obj <- obj
  u_obj$u_hat <- rev(obj$u_hat)
  expect_error(efficiency_ci(u_obj), "u_hat")

  ## Wrong length is still refused, and the message names the count.
  short <- obj
  short$exp_u_hat <- obj$exp_u_hat[1:3]
  expect_error(efficiency_ci(short), "cannot be matched")
  expect_error(efficiency_ci(short), "3 values")
})

test_that("a stored predictor of the wrong length is refused by name", {
  skip_on_cran()
  fit <- a69_fit()
  n <- length(fit$u_posterior$mu_star)

  ## Half length divides n, which is exactly the case a data frame recycles
  ## instead of rejecting (issue #2).
  short <- fit
  short$exp_u_hat <- as.numeric(fit$exp_u_hat)[seq_len(n / 2)]
  expect_error(efficiency_ci(short), "exp_u_hat")
  expect_error(efficiency_ci(short), "cannot be matched")
  expect_error(efficiency_ci(short), as.character(n / 2))

  short_u <- fit
  short_u$u_hat <- as.numeric(sfa:::.jlms_u(fit$u_posterior$mu_star,
                                           fit$u_posterior$sigma_star))[seq_len(n / 2)]
  expect_error(efficiency_ci(short_u), "u_hat")
  expect_error(efficiency_ci(short_u), "cannot be matched")
})

test_that("a stored predictor of the right length and wrong contents is refused", {
  skip_on_cran()
  fit <- a69_fit()
  stored <- as.numeric(fit$exp_u_hat)

  ## Reversing a fitted fixture is only a corruption if its scores actually
  ## differ across observations, which is a property of where the optimizer
  ## landed. Skip rather than fail if it collapsed -- the handmade block above
  ## already pins the guard's sensitivity with no such dependence.
  spread <- max(abs(rev(stored) - stored))
  skip_if_not(spread > 1e-2,
              sprintf("fitted scores too flat for reversal to corrupt (%.2g)", spread))

  ## Reversal is the sharp case: same length, same values, same mean, same
  ## range -- only the pairing with observations is destroyed. Nothing
  ## downstream can see it.
  rev_fit <- fit
  rev_fit$exp_u_hat <- rev(stored)
  expect_error(efficiency_ci(rev_fit), "exp_u_hat")
  expect_error(efficiency_ci(rev_fit), "do not belong to these observations")

  ## A single perturbed observation, far smaller than a reversal, is still
  ## caught -- the guard is an identity, not a tolerance on aggregates.
  one_off <- fit
  v <- stored
  v[7] <- v[7] + 1e-4
  one_off$exp_u_hat <- v
  expect_error(efficiency_ci(one_off), "do not belong to these observations")

  ## u_hat gets the same treatment when the fit carries it (NHN does not, so
  ## set it to the correct value first and then corrupt it).
  u_fit <- fit
  u_fit$u_hat <- as.numeric(sfa:::.jlms_u(fit$u_posterior$mu_star,
                                          fit$u_posterior$sigma_star))
  expect_s3_class(efficiency_ci(u_fit), "data.frame")
  u_fit$u_hat <- rev(u_fit$u_hat)
  expect_error(efficiency_ci(u_fit), "u_hat")
})

test_that("the guard does not fire on fits the package itself produced", {
  skip_on_cran()
  ## The guard is only worth having if it is silent on real fits. Both
  ## orientations, every model with a truncated-normal posterior.
  d <- cs_small(N = 150)
  specs <- list(
    list(nm = "NHN", f = y_pcs ~ x1 + x2),
    list(nm = "NHN_Z", f = y_pcs_z ~ x1 + x2 | z),
    list(nm = "NE", f = y_pcs ~ x1 + x2),
    list(nm = "NTN", f = y_pcs_tn ~ x1 + x2)
  )
  for (sp in specs) {
    for (ind in c(TRUE, FALSE)) {
      fit <- suppressWarnings(
        sfm(sp$f, model_name = sp$nm, data = d, inefdec = ind)
      )
      lab <- paste(sp$nm, ind)
      expect_s3_class(efficiency_ci(fit), "data.frame")
      ## And the stored field is genuinely there to be checked, so the pass is
      ## not a vacuous skip over a NULL.
      expect_false(is.null(fit$exp_u_hat), info = lab)
    }
  }
})

test_that("a fit with the field removed is still accepted", {
  skip_on_cran()
  ## The guard checks what is present; it must not turn an absent field into an
  ## error, which is the population A43's own guard gap depends on.
  fit <- a69_fit()
  stripped <- fit
  stripped$exp_u_hat <- NULL
  expect_s3_class(efficiency_ci(stripped), "data.frame")
  expect_equal(nrow(efficiency_ci(stripped)),
               length(fit$u_posterior$mu_star))
})
