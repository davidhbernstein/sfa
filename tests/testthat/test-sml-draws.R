## Draws for simulated maximum likelihood, .sml_draws() and sfm()'s options.
##
## The properties pinned here are the ones Train (2002, ch. 9) requires and the
## ones that were previously established inline in sfm.R and would be easy to
## lose in a refactor:
##   - each unit gets its own CONTIGUOUS block, not a strided subsequence;
##   - leading elements are discarded;
##   - draws never touch 0 or 1;
##   - antithetics are exact mirrors on the uniform scale;
##   - and the default reproduces the pre-existing construction EXACTLY, so no
##     fitted result moves.

test_that("the default reproduces the previous construction byte for byte", {
  for (cfg in list(c(500, 100), c(377, 213), c(120, 64))) {
    n <- cfg[1]; Ns <- cfg[2]
    hseq <- randtoolbox::halton(n * Ns + 1000, 1, start = 1, normal = FALSE)[-c(1:1000)]
    old <- matrix(pmin(pmax(hseq, 1e-6), 1 - 1e-6), nrow = n, ncol = Ns, byrow = TRUE)
    new <- sfa:::.sml_draws(n, Ns, dim = 1L, sim_type = "halton",
                            antithetics = FALSE, burn = 1000, clamp = 1e-6)[[1]]
    expect_identical(old, new, info = paste(n, Ns))
  }
})

test_that("every unit gets a contiguous block, not a strided subsequence", {
  ## The failure this guards against: a column-major fill gives unit 1 the
  ## elements h_1, h_{1+n}, h_{1+2n}, ... of a van der Corput sequence, which
  ## for n = 500, n_draws = 100 spans only [0.50, 0.75]. A contiguous block
  ## spans essentially the whole interval.
  M <- sfa:::.sml_draws(500, 100, dim = 1L, sim_type = "halton")[[1]]
  expect_equal(dim(M), c(500L, 100L))
  spans <- apply(M, 1, function(r) diff(range(r)))
  expect_gt(min(spans), 0.9)
  ## And rows must differ from one another -- sharing one block across units is
  ## exactly the defect this constructor exists to prevent.
  expect_false(isTRUE(all.equal(M[1, ], M[2, ])))
  expect_false(isTRUE(all.equal(M[1, ], M[nrow(M), ])))
})

test_that("burn-in is applied and draws stay strictly inside (0, 1)", {
  no_burn <- sfa:::.sml_draws(10, 20, sim_type = "halton", burn = 0)[[1]]
  burned  <- sfa:::.sml_draws(10, 20, sim_type = "halton", burn = 1000)[[1]]
  expect_false(isTRUE(all.equal(no_burn, burned)))
  for (M in list(no_burn, burned)) {
    expect_true(all(M > 0 & M < 1))
    expect_true(all(is.finite(qnorm(M))))
  }
})

test_that("antithetics are exact mirrors and preserve the requested count", {
  M <- sfa:::.sml_draws(8, 20, sim_type = "halton", antithetics = TRUE)[[1]]
  expect_equal(dim(M), c(8L, 20L))
  ## First half and second half must sum to 1 elementwise.
  expect_equal(M[, 1:10] + M[, 11:20], matrix(1, 8, 10), tolerance = 1e-12)
  ## A symmetric transform then gives exact sign mirrors, which is the point:
  ## qnorm(u) = -qnorm(1 - u).
  Z <- qnorm(M)
  expect_equal(Z[, 1:10], -Z[, 11:20], tolerance = 1e-10)
  ## An odd count must still return exactly what was asked for.
  Modd <- sfa:::.sml_draws(4, 7, sim_type = "halton", antithetics = TRUE)[[1]]
  expect_equal(ncol(Modd), 7L)
})

test_that("each sequence type runs and gives distinct, valid draws", {
  got <- lapply(c("halton", "sobol", "torus", "uniform"), function(st)
    sfa:::.sml_draws(30, 40, sim_type = st, seed = 7)[[1]])
  names(got) <- c("halton", "sobol", "torus", "uniform")
  for (nm in names(got)) {
    expect_equal(dim(got[[nm]]), c(30L, 40L), info = nm)
    expect_true(all(got[[nm]] > 0 & got[[nm]] < 1), info = nm)
  }
  ## They must actually be different sequences.
  expect_false(isTRUE(all.equal(got$halton, got$sobol)))
  expect_false(isTRUE(all.equal(got$halton, got$torus)))
  expect_false(isTRUE(all.equal(got$halton, got$uniform)))
})

test_that("multiple dimensions come back as separate, uncorrelated-enough blocks", {
  ## Portability contract for psfm(), which integrates in 2 dimensions: the
  ## constructor must return one matrix per dimension, each blocked by unit.
  D <- sfa:::.sml_draws(50, 60, dim = 2L, sim_type = "halton")
  expect_length(D, 2L)
  for (M in D) expect_equal(dim(M), c(50L, 60L))
  expect_false(isTRUE(all.equal(D[[1]], D[[2]])))
  ## Primes 2 and 3 with a burn-in leave the dimensions near-uncorrelated
  ## already, which is why the old brute-force decorrelation was unnecessary.
  expect_lt(abs(cor(as.vector(D[[1]]), as.vector(D[[2]]))), 0.05)
})

test_that("a seed makes the draws reproducible and restores the RNG stream", {
  set.seed(123); before <- runif(1)
  set.seed(123)
  a <- sfa:::.sml_draws(10, 10, sim_type = "uniform", seed = 42)[[1]]
  after <- runif(1)
  b <- sfa:::.sml_draws(10, 10, sim_type = "uniform", seed = 42)[[1]]
  expect_equal(a, b)
  ## The caller's stream must be where it would have been.
  expect_equal(after, before)
})

test_that("sfm exposes the options and its default is unchanged", {
  skip_on_cran()
  d <- cs_small(N = 250)
  base <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d))
  same <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d,
                               sim_type = "halton", antithetics = FALSE))
  expect_equal(coef(base), coef(same), tolerance = 1e-10)

  ## A different scheme must actually change the fit, or the argument is inert.
  alt <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d,
                              sim_type = "sobol", sim_scrambling = 1L, sim_seed = 5))
  expect_false(isTRUE(all.equal(unname(coef(base)), unname(coef(alt)))))
  expect_true(all(is.finite(coef(alt))))
})

test_that("bad simulation arguments are rejected", {
  expect_error(sfa:::.sml_draws(0, 10), "at least 1")
  ## A NULL or NA size must say so, not fall through to a comparison error.
  expect_error(sfa:::.sml_draws(10, 10, dim = NULL), "single positive")
  expect_error(sfa:::.sml_draws(NA, 10), "single positive")
  expect_error(sfa:::.sml_draws(10, 10, sim_type = "sausage"), "should be one of")
  skip_on_cran()
  d <- cs_small(N = 120)
  expect_error(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d, antithetics = "yes"),
               "TRUE or FALSE")
  expect_error(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d, sim_burn = -1),
               "non-negative")
})
