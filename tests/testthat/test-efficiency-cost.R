## Gap A41, from GitHub issue #2. efficiency(logDepVar = FALSE) formed the
## production ratio 1 - u/f for every fit, including cost frontiers, where
## y = f + u + v and efficiency is minimum over actual cost, f / (f + u).

.lvl_data <- function() {
  set.seed(11)
  n <- 400
  x <- runif(n, 1, 3)
  v <- rnorm(n, 0, 0.3)
  u <- abs(rnorm(n, 0, 0.6))
  data.frame(x = x, y_prod = 5 + 2 * x + v - u, y_cost = 5 + 2 * x + v + u)
}

test_that("level-scale efficiency of a cost frontier is f / (f + u), in (0, 1]", {
  d <- .lvl_data()
  fc <- suppressWarnings(sfm(y_cost ~ x, model_name = "NHN", data = d, inefdec = FALSE))
  f <- fitted(fc)
  u <- -log(efficiency(fc, type = "jlms"))
  te <- efficiency(fc, type = "jlms", logDepVar = FALSE)
  expect_equal(te, f / (f + u))
  expect_true(all(te > 0 & te <= 1))
})

test_that("the production frontier is unchanged: 1 - u / f", {
  d <- .lvl_data()
  fp <- suppressWarnings(sfm(y_prod ~ x, model_name = "NHN", data = d))
  f <- fitted(fp)
  u <- -log(efficiency(fp, type = "jlms"))
  expect_equal(efficiency(fp, type = "jlms", logDepVar = FALSE), 1 - u / f)
})

test_that("inefdec passed as a variable is read, and an unreadable one is refused", {
  d <- .lvl_data()
  flag <- FALSE
  fc <- suppressWarnings(sfm(y_cost ~ x, model_name = "NHN", data = d, inefdec = FALSE))
  fv <- suppressWarnings(sfm(y_cost ~ x, model_name = "NHN", data = d, inefdec = flag))
  expect_equal(efficiency(fv, logDepVar = FALSE), efficiency(fc, logDepVar = FALSE))
  expect_equal(predict(fv, type = "response"), predict(fc, type = "response"))
  ## A cost frontier's prediction adds inefficiency to the frontier. NHN fits
  ## store exp_u_hat, not u_hat, so predict() uses -log(exp_u_hat).
  u <- -log(pmax(fc$exp_u_hat, .Machine$double.xmin))
  expect_equal(as.numeric(predict(fc, type = "response")), as.numeric(fitted(fc) + u))

  fx <- fv
  fx$call$inefdec <- as.name("no_such_variable_anywhere")
  expect_error(efficiency(fx, logDepVar = FALSE), "cannot tell whether")
  expect_error(predict(fx, type = "response"), "cannot tell whether")
})

test_that("efficiency_ci() refuses point predictors of the wrong length", {
  ## A data frame column silently recycles a vector whose length divides the row
  ## count, which would pair estimates with the wrong observations (issue #2).
  d <- .lvl_data()
  f <- suppressWarnings(sfm(y_prod ~ x, model_name = "NHN", data = d))
  expect_silent(efficiency_ci(f))
  g <- f
  g$u_hat <- rep(0.1, nrow(d) / 2)
  expect_error(efficiency_ci(g), "`object\\$u_hat` has 200 values")
  h <- f
  h$exp_u_hat <- h$exp_u_hat[seq_len(nrow(d) / 2)]
  expect_error(efficiency_ci(h), "`object\\$exp_u_hat` has 200 values")
})
