## Gap A46: psfm_bootstrap() reproducibility.
##
## Replication b is seeded inside the worker, so the draws cannot depend on
## how many workers there are or on the caller's RNG. What was wrong was
## seed_offset: the seed was b + seed_offset, so a run at offset 1 repeated
## replications 2..BOOT of the run at offset 0 exactly -- measured, 5 of 6 --
## although ?psfm_bootstrap promised "reproducible-but-distinct seeds across
## multiple bootstrap runs". Two runs pooled as independent shared nearly every
## draw, and the pooled standard error understated the spread.

.a46_fit <- function() {
  d <- as.data.frame(data_gen_p(t = 4, N = 30, rand = 11, sig_u = 1,
    sig_v = 0.3, sig_r = 0.4, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  suppressWarnings(psfm(y_tfe ~ x1 + x2, "TFE_WMLE", d, individual = "name"))
}
.a46_boot <- function(f, cores, ...) {
  suppressWarnings(psfm_bootstrap(f, numCores = cores, BOOT = 4,
    individual = "name", ...))
}
.a46_par <- function(b) unname(b$boot_par[, colnames(b$boot_par) != "hours"])

test_that("runs at different seed_offset share no replication", {
  skip_on_cran()
  f <- .a46_fit()
  p0 <- .a46_par(.a46_boot(f, 1, seed_offset = 0))
  p1 <- .a46_par(.a46_boot(f, 1, seed_offset = 1))
  expect_true(all(is.finite(p0)))
  shared <- outer(seq_len(nrow(p0)), seq_len(nrow(p1)),
    Vectorize(function(i, j) isTRUE(all.equal(p0[i, ], p1[j, ], tolerance = 1e-12))))
  expect_false(any(shared))
})

test_that("the draws do not depend on the worker count or the caller's RNG", {
  skip_on_cran()
  f <- .a46_fit()
  set.seed(99)
  s0 <- .Random.seed
  b1 <- .a46_boot(f, 1)
  expect_identical(.Random.seed, s0)
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  set.seed(5)
  b2 <- .a46_boot(f, 2)
  expect_identical(.a46_par(b1), .a46_par(b2))
  expect_identical(b1$boot_eff, b2$boot_eff)
})

test_that("seed_offset is validated", {
  f <- structure(list(out = matrix(0), data = data.frame(name = 1),
    formula = y ~ x, model_name = "TRE"), class = "sfareg")
  expect_error(psfm_bootstrap(f, 1, 2, "name", seed_offset = 0.5, inefdec = TRUE),
    "single whole number")
  expect_error(psfm_bootstrap(f, 1, 100001, "name", inefdec = TRUE), "exceeds")
  expect_error(psfm_bootstrap(f, 1, 2, "name", seed_offset = 1e5, inefdec = TRUE),
    "must lie between")
})
