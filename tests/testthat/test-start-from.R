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

test_that("seeding never makes the fit worse, and sometimes rescues it", {
  skip_on_cran()
  ## Seed 4 is the case the feature exists for: from the default start NG lands
  ## on an optimum 176 log-likelihood points below the one reachable from the
  ## NE fit. Most samples are unaffected, which is why this is plumbing rather
  ## than a change of estimator -- but "never worse" is the property that makes
  ## it safe to reach for.
  d <- sf_data(4)
  ne <- sfm(y ~ x1 + x2, model_name = "NE", data = d)
  a <- sfm(y ~ x1 + x2, model_name = "NG", data = d)
  b <- suppressMessages(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = ne)
  )
  expect_gt(-b$opt$value, -a$opt$value + 50)
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
