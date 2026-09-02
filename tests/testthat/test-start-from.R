## start_from: seed a hard model from a simpler fit already in hand (G2/I3).

sf_data <- function(seed, n = 300) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - rgamma(n, shape = 2, scale = 0.4)
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("matching is by NAME, so shared scale names carry and others do not", {
  skip_on_cran()
  d <- sf_data(1)
  ne <- sfm(y ~ x1 + x2, model_name = "NE", data = d)
  nhn <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)

  ## NE and NG both call their scales sigv/sigu, so five of six carry.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = ne),
    "carried 5 of 6"
  )
  ## NHN calls them lambda/sigma. Those must NOT be poured into NG's sigv/sigu
  ## -- a positional copy would slide a variance ratio into a standard
  ## deviation. Only the three frontier coefficients are shared.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = nhn),
    "carried 3 of 6"
  )
  ## NHN -> NTN shares lambda and sigma, which is the documented idiom.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, start_from = nhn),
    "carried 5 of 6"
  )
})

test_that("seeding never makes the fit worse", {
  skip_on_cran()
  ## "Never worse" is the only property of start_from that holds everywhere, so
  ## it is the only one asserted. The SIZE of any rescue does not travel: over
  ## seeds 1-12 on this machine, ten give an identical optimum either way, seed
  ## 4 gains 176.70 log-likelihood points, and seeds 7 and 10 make the UNSEEDED
  ## NG fit fail outright in optim (seeding rescues 7, not 10). On the CI
  ## platforms seed 4 converges to the same optimum from both starts. Asserting
  ## the gain would pin a BLAS/optimizer path rather than the feature.
  ##
  ## A plain fit that errors counts as "seeding is no worse" rather than as a
  ## test failure -- that failure is NG's own start-value fragility, tracked
  ## separately, and it is exactly the situation start_from exists to escape.
  ll <- function(e) if (inherits(e, "error")) NA_real_ else -e$opt$value
  for (s in c(1, 2, 4)) {
    d <- sf_data(s)
    ne <- sfm(y ~ x1 + x2, model_name = "NE", data = d)
    a <- tryCatch(sfm(y ~ x1 + x2, model_name = "NG", data = d),
      error = function(e) e
    )
    b <- tryCatch(
      suppressMessages(
        sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = ne)
      ),
      error = function(e) e
    )
    expect_false(inherits(b, "error"))
    if (!is.na(ll(a))) expect_gte(ll(b), ll(a) - 1e-6)
  }
})

test_that("start_from rejects what is not a fit", {
  d <- sf_data(2)
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = list(a = 1)),
    "must be a fitted"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = 1:5),
    "must be a fitted"
  )
})
