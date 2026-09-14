## Draw controls for psfm()'s simulated-ML panel models (gap H10). The draws
## come from .sml_draws(), shared with sfm(); these tests pin the panel wrapper
## and the arguments psfm() passes through to it.

test_that(".gtre_halton_draws() returns one finite R x 2 block per firm for every sim_type", {
  for (st in c("halton", "sobol", "torus", "uniform")) {
    dl <- .gtre_halton_draws(N = 7, R = 20, rand.gtre = 3, sim_type = st)
    expect_length(dl, 7)
    expect_true(all(vapply(dl, function(m) identical(dim(m), c(20L, 2L)), logical(1))), info = st)
    expect_true(all(is.finite(unlist(dl))), info = st)
  }
})

test_that(".gtre_halton_draws() is reproducible, leaves the RNG alone, and sim_type matters", {
  set.seed(5)
  before <- .Random.seed
  a <- .gtre_halton_draws(N = 5, R = 30, rand.gtre = 11, sim_type = "sobol", scrambling = 0L)
  expect_identical(.Random.seed, before)
  b <- .gtre_halton_draws(N = 5, R = 30, rand.gtre = 11, sim_type = "sobol", scrambling = 0L)
  expect_identical(a, b)
  h <- .gtre_halton_draws(N = 5, R = 30, rand.gtre = 11)
  expect_false(isTRUE(all.equal(a, h)))
  anti <- .gtre_halton_draws(N = 5, R = 30, rand.gtre = 11, antithetics = TRUE)
  expect_false(isTRUE(all.equal(anti, h)))
})

test_that("psfm() validates the draw controls before fitting", {
  d <- panel_small(t = 4, N = 20)
  expect_error(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d, individual = "name",
                    antithetics = NA), "antithetics")
  expect_error(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d, individual = "name",
                    sim_burn = -1), "sim_burn")
  expect_error(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d, individual = "name",
                    sim_type = "ghalton"))
})

test_that("psfm() passes the draw controls to the simulated-ML fits", {
  skip_on_cran()
  d <- panel_small(t = 5, N = 30)
  fit <- function(...) {
    suppressWarnings(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d,
                          individual = "name", halton_num = 30, rand.gtre = 4,
                          maxit.bobyqa = 100, maxit.optim = 50, ...))
  }
  default <- fit()
  explicit <- fit(sim_type = "halton", antithetics = FALSE, sim_burn = 1000)
  expect_identical(default$out, explicit$out)

  ## Sobol's first dimension is Halton's (both are base-2 van der Corput), so
  ## on this design, where sigma_r sits at its bound and the second dimension
  ## barely enters, a Sobol fit differs from the Halton one only around 1e-9.
  ## Uniform draws are the check that the choice reaches the likelihood.
  sob1 <- fit(sim_type = "sobol")
  sob2 <- fit(sim_type = "sobol")
  expect_true(all(is.finite(sob1$coefficients)))
  expect_identical(sob1$out, sob2$out)
  expect_false(identical(sob1$opt$value, default$opt$value))
  unif <- fit(sim_type = "uniform")
  expect_false(isTRUE(all.equal(unif$opt$value, default$opt$value)))
})
