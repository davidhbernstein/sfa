## psfm("LS")'s convergence flag (gap A83).
##
## `.fit_ls93()` iterates alternating least squares and BREAKS when its
## tolerance is met, then reported `converged <- it < maxit`. A run that met
## its tolerance ON the last allowed iteration breaks with `it == maxit`, so it
## was reported as NOT converged and warned the user to treat `delta` and the
## inefficiency scores as provisional -- on a fit that had in fact converged.
##
## Asserted by calling `.fit_ls93()` directly, for two reasons. `psfm()` never
## passes `maxit` through to it (R/psfm.R calls `.fit_ls93(setup, inefdec)`
## with the default 500), so the boundary is not reachable through the public
## interface at all; and driving `maxit` to either side of the iteration count
## this design actually needs tests the FLAG rather than where the iteration
## happens to land.

ls_setup <- function(seed = 4, N = 30, t = 5) {
  d <- data_gen_p(N = N, t = t, rand = seed, sig_u = 1, sig_v = .5,
                  sig_r = .5, sig_h = .5, cons = .5, beta1 = .5, beta2 = .5)
  .classical_panel_setup(y_ssfe ~ x1 + x2, d, "name", NULL, "y_ssfe", "LS")
}

test_that("a converged LS fit reports converged and does not warn", {
  skip_on_cran()
  s <- ls_setup()
  expect_silent(f <- .fit_ls93(s, inefdec = TRUE))
  expect_true(isTRUE(f$converged))
  expect_lt(f$iterations, 500L)
})

test_that("the flag tracks the tolerance, not the iteration count", {
  skip_on_cran()
  s <- ls_setup()
  need <- .fit_ls93(s, inefdec = TRUE)$iterations
  skip_if(need < 2L, "this design converges too fast to place the boundary")

  ## maxit set to EXACTLY the number of iterations needed, so the break lands
  ## on the final allowed pass. This is the case the old flag got wrong: it
  ## reported converged = FALSE and warned.
  expect_silent(g <- .fit_ls93(s, inefdec = TRUE, maxit = need))
  expect_equal(g$iterations, need)
  expect_true(isTRUE(g$converged))

  ## One fewer and it genuinely cannot converge, so the flag must flip and the
  ## warning must fire. A flag hard-coded TRUE would pass everything above and
  ## fail here.
  expect_warning(h <- .fit_ls93(s, inefdec = TRUE, maxit = need - 1L),
                 "without meeting its tolerance")
  expect_false(isTRUE(h$converged))
})

test_that("the flag reaches the fitted object, under $ls_converged", {
  skip_on_cran()
  d <- data_gen_p(N = 30, t = 5, rand = 4, sig_u = 1, sig_v = .5, sig_r = .5,
                  sig_h = .5, cons = .5, beta1 = .5, beta2 = .5)
  f <- psfm(y_ssfe ~ x1 + x2, data = d, model_name = "LS", individual = "name")
  expect_s3_class(f, "sfareg")
  ## `$ls_converged`, not `$converged` (R/psfm.R surfaces it under that name);
  ## asserted so the next reader does not look for the wrong slot, as this test
  ## first did.
  expect_true(isTRUE(f$ls_converged))
  expect_null(f$converged)
})
