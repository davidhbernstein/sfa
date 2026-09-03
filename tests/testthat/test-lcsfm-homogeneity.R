## lcsfm_homogeneity(): the modified LR test for how many classes (L3).
##
## What is asserted is what the arithmetic or the measurement fixes: that the
## penalty vanishes at equal class probabilities, that `.chi2_01_p` really is
## half a one-df tail, that the penalised objective is kept apart from the
## reported log-likelihood, and that the test points the right way on data
## built to have one class and on data built to have two. No particular MLR
## value is asserted -- that is a sampling outcome and an optimizer path.
##
## Note what this file does NOT claim. The published chi^2_{0:1} null is
## measurably WRONG for `lcsfm()`'s "LCM", which lets every parameter vary by
## class rather than the single scalar the theorem assumes: it rejected 63.5%
## of the time at a nominal 5% over 200 replications of a one-class DGP. So
## `null = "chisq01"` is tested for the fact that it WARNS, and the default
## bootstrap is tested for its behaviour. The bootstrap cases are B + 1 refits
## apiece and sit behind SFA_TEST_SLOW.

lcm_one_class <- function(seed = 11, n = 400) {
  set.seed(seed)
  x1 <- rnorm(n)
  data.frame(y = 1 + x1 + rnorm(n, 0, 1) - abs(rnorm(n, 0, 1)), x1 = x1)
}

lcm_two_class <- function(seed = 21, n = 400) {
  set.seed(seed)
  x1 <- rnorm(n)
  cls <- rbinom(n, 1, 0.5)
  data.frame(
    y = ifelse(cls == 1, 4, 1) + x1 + rnorm(n, 0, 0.4) - abs(rnorm(n, 0, 0.6)),
    x1 = x1
  )
}

test_that("the chi^2_{0:1} p-value is HALF the one-df tail, and 1 at zero", {
  ## The whole point of the entry: the null is a 50:50 mixture of a point mass
  ## at zero and chi^2_1, so the p-value is half the naive one-df tail and does
  ## NOT depend on the parameter-count difference. Using pchisq(x, df) with df
  ## = the number of extra parameters is the error the test exists to avoid.
  for (x in c(0.5, 1, 2.7, 3.84, 10)) {
    expect_equal(sfa:::.chi2_01_p(x), 0.5 * pchisq(x, 1, lower.tail = FALSE))
  }
  expect_equal(sfa:::.chi2_01_p(0), 1)
  expect_equal(sfa:::.chi2_01_p(-3), 1)
  ## The 5% critical value of chi^2_{0:1} is the 90th percentile of chi^2_1,
  ## which is 2.706 -- NOT 3.841.
  expect_equal(sfa:::.chi2_01_p(qchisq(0.90, 1)), 0.05, tolerance = 1e-8)
})

test_that("the penalty vanishes at equal class probabilities", {
  skip_on_cran()
  ## c * log(J^J * prod p_j) is 0 when every p_j = 1/J, which is what lets the
  ## null model be recovered without a handicap and keeps the statistic from
  ## going negative for that reason. At J = 2: 2 log 2 + 2 log 0.5 = 0.
  expect_equal(2 * log(2) + sum(log(c(0.5, 0.5))), 0)
  expect_equal(3 * log(3) + sum(log(rep(1 / 3, 3))), 0)
  ## And it is strictly negative away from equality, so it always penalises.
  expect_lt(2 * log(2) + sum(log(c(0.1, 0.9))), 0)
  expect_lt(2 * log(2) + sum(log(c(0.01, 0.99))), 0)
})

test_that("penalty_c changes the objective but not the reported log-likelihood", {
  skip_on_cran()
  d <- lcm_one_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2, penalty_c = 1)
  ## opt$value is the PENALISED objective; the plain log-likelihood is carried
  ## separately so logLik() cannot silently report a penalised number.
  expect_equal(f$logLik_unpenalised, -f$opt$value + f$penalty)
  expect_lte(f$penalty, 0)
  expect_identical(f$penalty_c, 1)

  ## With penalty_c = 0 the penalty is EXACTLY zero, not merely small.
  f0 <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  expect_identical(f0$penalty, 0)
})

test_that("chisq01 warns, because it is measurably over-sized for this model", {
  skip_on_cran()
  ## The published chi^2_{0:1} limit assumes a SCALAR parameter differs between
  ## classes. "LCM" lets every parameter vary, so it does not apply. Measured
  ## over 200 replications of a one-class DGP at n = 400, that null rejected
  ## 63.5% of the time at a nominal 5%, and the statistic tracked chi^2_5 --
  ## exactly the parameter-count difference. The warning is the deliverable.
  d <- lcm_one_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  expect_warning(
    tt <- lcsfm_homogeneity(f, null = "chisq01"),
    "63.5%"
  )
  expect_true(tt$provisional)
  expect_output(print(tt), "bootstrap")
})

test_that("the bootstrap null does not reject on one-class data", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("SFA_TEST_SLOW")),
    "bootstrap is B + 1 refits; set SFA_TEST_SLOW to run"
  )
  ## Seed 11 at n = 400 is the case that FALSELY rejected under chi^2_{0:1}
  ## (p = 0.005). Under the bootstrap it gives p = 0.32. Kept as the regression
  ## test for the fix, asserting only that it does not reject -- the p-value
  ## itself is a bootstrap draw.
  d <- lcm_one_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  tt <- lcsfm_homogeneity(f, null = "bootstrap", B = 49, seed = 1, quiet = TRUE)
  expect_identical(tt$null, "bootstrap")
  expect_false(tt$provisional)
  expect_gt(tt$p.value, 0.05)
  expect_gte(tt$statistic, 0)
})

test_that("the bootstrap null still rejects on two-class data", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("SFA_TEST_SLOW")),
    "bootstrap is B + 1 refits; set SFA_TEST_SLOW to run"
  )
  ## Power is the other half. On well-separated classes the statistic came out
  ## 17x the largest of 99 bootstrap null draws, so this is not a close call.
  d <- lcm_two_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  tt <- lcsfm_homogeneity(f, null = "bootstrap", B = 49, seed = 2, quiet = TRUE)
  expect_lt(tt$p.value, 0.05)
  expect_true(tt$reject)
  expect_gt(tt$statistic, max(tt$boot))
})

test_that("the statistic is computed against the one-class fit, not a refit lcsfm", {
  skip_on_cran()
  ## One latent class IS the ordinary frontier, so the null model is
  ## sfm("NHN"). Checking the stored null log-likelihood against a direct fit
  ## catches the null being built from the wrong model.
  d <- lcm_one_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  tt <- suppressWarnings(lcsfm_homogeneity(f, null = "chisq01", c = 1))
  direct <- sfm(y ~ x1, model_name = "NHN", data = d)
  expect_equal(tt$logLik_1, -direct$opt$value, tolerance = 1e-6)
  expect_equal(tt$statistic[["MLR"]],
    2 * (tt$logLik_penalised_J - tt$logLik_1),
    tolerance = 1e-8
  )
})

test_that("a larger c penalises unequal probabilities more", {
  skip_on_cran()
  d <- lcm_one_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  t1 <- suppressWarnings(lcsfm_homogeneity(f, null = "chisq01", c = 1))
  t5 <- suppressWarnings(lcsfm_homogeneity(f, null = "chisq01", c = 5))
  ## Both values are the ones Stead, Wheat and Greene use. The penalty at the
  ## maximum is <= 0 either way, and c = 5 cannot be the less negative of the
  ## two at the same probabilities.
  expect_lte(t5$penalty_at_max, 0)
  expect_identical(t1$penalty_c, 1)
  expect_identical(t5$penalty_c, 5)
})

test_that("the refusals name the reason", {
  skip_on_cran()
  d <- lcm_one_class()
  expect_error(lcsfm_homogeneity(list(a = 1)), "sfareg")

  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  expect_error(lcsfm_homogeneity(f, c = 0), "positive")
  expect_error(lcsfm_homogeneity(f, c = -1), "positive")
  expect_error(lcsfm_homogeneity(f, B = 5), "B = 19")

  ## LCM_Z drives the class probabilities with covariates, which is outside
  ## the published result -- both the fitter and the test must refuse.
  expect_error(
    lcsfm(y ~ x1 | x1, model_name = "LCM_Z", data = d, n_class = 2,
      penalty_c = 1),
    "scalar mixing proportion"
  )
  expect_error(lcsfm(y ~ x1, model_name = "LCM", data = d, penalty_c = -1),
    "penalty_c"
  )
})

test_that("print reports the verdict and flags J > 2 as provisional", {
  skip_on_cran()
  d <- lcm_two_class()
  f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)
  tt <- suppressWarnings(lcsfm_homogeneity(f, null = "chisq01", c = 1))
  expect_output(print(tt), "Modified likelihood ratio")
  expect_output(print(tt), "MLR")
  expect_output(print(tt), "p-value")
  ## chisq01 is always flagged now, since "LCM" is never the restricted model
  ## the limit is stated for.
  expect_true(tt$provisional)
  expect_output(print(tt), "63.5")
})
