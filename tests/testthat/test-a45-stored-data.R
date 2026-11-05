## Gap A45: a fit carries its data, so the methods work once the data is out
## of scope -- the case no environment search can reach is a formula built in
## one place and a fit made in another, on data local to the fitting function.

a45_data <- function() {
  set.seed(45)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 5 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.3) - abs(rnorm(n, sd = 0.6))
  d <- data.frame(y, x1, x2)
  d$y[c(3, 7)] <- NA
  d
}

test_that("methods work on a fit whose data is local to another function (A45)", {
  frm <- y ~ x1 + x2
  fit_away <- function(fitter) {
    local_only_a45 <- a45_data()
    suppressWarnings(fitter(frm, local_only_a45))
  }
  fits <- list(
    sfm = function(f, d) sfm(f, model_name = "NHN", data = d),
    zsfm = function(f, d) zsfm(f, model_name = "ZISF", data = d)
  )
  d <- a45_data()
  for (nm in names(fits)) {
    away <- fit_away(fits[[nm]])
    expect_false(exists("local_only_a45"))
    expect_identical(away$data, d)
    here <- suppressWarnings(fits[[nm]](frm, d))
    expect_equal(fitted(away), fitted(here), info = nm)
    expect_equal(residuals(away), residuals(here), info = nm)
    expect_equal(predict(away), predict(here), info = nm)
    expect_identical(nobs(away), 198L, info = nm)
  }
  away <- fit_away(fits$sfm)
  here <- suppressWarnings(fits$sfm(frm, d))
  expect_equal(efficiency(away, logDepVar = FALSE), efficiency(here, logDepVar = FALSE))
})

test_that(".sfa_keep_data() never replaces data a branch stored itself", {
  own <- data.frame(a = 1:3)
  r <- .sfa_keep_data(list(data = own), data.frame(b = 1:5))
  expect_identical(r$data, own)
  ## an exact match: a component that merely starts with "data" does not count
  r <- .sfa_keep_data(list(data_i = 1), own)
  expect_identical(r[["data", exact = TRUE]], own)
  expect_null(.sfa_keep_data(list(), NULL)[["data", exact = TRUE]])
})
