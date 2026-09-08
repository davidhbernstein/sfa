## Endogenous regressors: ivsfm(), Amsler, Prokhorov and Schmidt (2016, 2017).

## APS (2016) Eq (10) and its reduced form. The reduced-form error eta is
## correlated with the NOISE v, not with the inefficiency u -- that is the
## model's defining assumption, so the generator has to respect it.
iv_data <- function(seed, n = 800, rho = 0.6, su = 1.0, sv = 0.5,
                    b = c(0.5, 0.8, -0.6), cost = FALSE, with_q = FALSE) {
  set.seed(seed)
  w1 <- rnorm(n)
  w2 <- rnorm(n)
  x1 <- rnorm(n)
  eta <- rnorm(n)
  v <- sv * (rho * eta + sqrt(1 - rho^2) * rnorm(n))
  x2 <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
  q <- if (with_q) rnorm(n) else NULL
  su_i <- if (with_q) su * exp(0.3 * q) else su
  u <- abs(rnorm(n, 0, su_i))
  y <- b[1] + b[2] * x1 + b[3] * x2 + v + if (cost) u else -u
  d <- data.frame(y = y, x1 = x1, x2 = x2, w1 = w1, w2 = w2)
  if (with_q) d$q <- q
  d
}

test_that("ivsfm returns a well-formed sfareg for each estimator", {
  d <- iv_data(1)
  for (mn in c("IVLIML", "IVCF", "C2SLS")) {
    f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = d, model_name = mn
    )
    expect_s3_class(f, "sfareg")
    expect_identical(f$model_name, mn)

    ## p x 3, one ROW per parameter, as everywhere else in the package.
    expect_equal(ncol(f$out), 3L)
    expect_identical(colnames(f$out), c("par", "st_err", "t-val"))
    expect_equal(unname(f$coefficients), unname(f$out[, "par"]))

    expect_true(all(c("sigma_u", "sigma_v") %in% rownames(f$out)))
    expect_gt(f$out["sigma_u", "par"], 0)
    expect_gt(f$out["sigma_v", "par"], 0)

    expect_equal(length(f$jlms), nrow(d))
    expect_true(all(f$jlms >= 0))
    expect_equal(f$efficiency, exp(-f$jlms))
    expect_equal(nobs(f), nrow(d))
  }
})

test_that("the ML estimators report rho and keep it inside the unit interval", {
  d <- iv_data(2)
  for (mn in c("IVLIML", "IVCF")) {
    f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = d, model_name = mn
    )
    expect_true("rho_x2" %in% rownames(f$out))
    ## sigma_c^2 = sigma_v^2 (1 - rho^2) has to stay positive, which the
    ## t/sqrt(1+t^2) parameterization makes unreachable rather than merely
    ## penalized. If this ever fails the transform has been broken.
    expect_lt(abs(f$rho), 1)
    expect_gt(f$sigma_c, 0)
    expect_lt(f$sigma_c, f$out["sigma_v", "par"] + 1e-8)
  }
  ## C2SLS never models the correlation it corrects for, so it reports none.
  fc <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "C2SLS"
  )
  expect_false("rho_x2" %in% rownames(fc$out))
})

test_that("ivsfm removes the endogeneity bias an uncorrected fit carries", {
  d <- iv_data(7, n = 2000)

  naive <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  b_naive <- naive$out[, "par"][["x2"]]

  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "IVLIML"
  )
  b_iv <- f$out["x2", "par"]

  ## The whole point of the entry point. Consistency alone would not show this:
  ## the uncorrected estimator is consistent for the WRONG number.
  expect_lt(abs(b_iv - (-0.6)), abs(b_naive - (-0.6)))
  expect_equal(unname(b_iv), -0.6, tolerance = 0.1)
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.1)
  expect_equal(unname(f$out["sigma_u", "par"]), 1.0, tolerance = 0.2)
  expect_equal(unname(f$out["rho_x2", "par"]), 0.6, tolerance = 0.2)
})

test_that("environmental variables give the 2017 model and recover delta", {
  d <- iv_data(3, n = 2000, with_q = TRUE)
  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    uhet = ~q, data = d, model_name = "IVLIML"
  )
  expect_true("delta_q" %in% rownames(f$out))
  ## sigma_u,i = sigma_u * exp(q_i * delta), delta = 0.3 in the generator.
  expect_equal(unname(f$out["delta_q", "par"]), 0.3, tolerance = 0.15)
  expect_equal(unname(f$out["x2", "par"]), -0.6, tolerance = 0.12)

  ## Without uhet the same data is fitted by the 2016 model, which has no delta.
  g <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "IVLIML"
  )
  expect_false(any(grepl("^delta", rownames(g$out))))
  expect_lt(nrow(g$out), nrow(f$out))
})

test_that("a cost frontier flips the sign of the inefficiency term", {
  d <- iv_data(4, n = 1500, cost = TRUE)
  fc <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "IVLIML", inefdec = FALSE
  )
  expect_equal(fc$S, -1)
  expect_equal(unname(fc$out["sigma_u", "par"]), 1.0, tolerance = 0.25)
  expect_equal(unname(fc$out["x2", "par"]), -0.6, tolerance = 0.12)
})

test_that("C2SLS is moment-based, so it reports no likelihood", {
  d <- iv_data(5)
  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "C2SLS"
  )
  ## Same footing as GTRE_SEQ1/GTRE_SEQ2: no $opt, so logLik() is NA with a
  ## warning rather than a number that would look like a comparable fit.
  expect_null(f$opt)
  expect_warning(ll <- logLik(f))
  expect_true(is.na(as.numeric(ll)))
})

test_that("the 2SLS residual uses the actual endogenous regressor", {
  ## APS (2016) Eq (11) is explicit that e = y - X b, NOT y - Xhat b. Using
  ## fitted values gives different second and third moments and so a different
  ## sigma_u. Verified by reconstructing both and checking which one the fit
  ## reproduces.
  d <- iv_data(6, n = 1200)
  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "C2SLS"
  )
  X <- model.matrix(~ x1 + x2, d)
  Z <- model.matrix(~ x1 + w1 + w2, d)
  b <- f$b_2sls
  Xh <- Z %*% solve(crossprod(Z), crossprod(Z, X))

  k3 <- sqrt(2 / pi) * (1 - 4 / pi)
  su_of <- function(e) (mean(e^3) / k3)^(1 / 3)
  su_right <- su_of(as.numeric(d$y - X %*% b))
  su_wrong <- su_of(as.numeric(d$y - Xh %*% b))

  expect_equal(unname(f$out["sigma_u", "par"]), su_right, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(su_right, su_wrong, tolerance = 1e-3)))
})

test_that("ivsfm rejects malformed calls", {
  d <- iv_data(8, n = 400)

  expect_error(ivsfm("y ~ x1", ~x2, ~ w1 + w2, d), "must be a formula")
  expect_error(ivsfm(~ x1 + x2, ~x2, ~ w1 + w2, d), "must have a left-hand side")
  expect_error(
    ivsfm(y ~ x1 + x2, ~x2, ~ w1 + w2, as.matrix(d)),
    "must be a data.frame"
  )
  expect_error(
    ivsfm(y ~ x1 | x2, ~x2, ~ w1 + w2, d),
    "must not contain a `\\|` segment"
  )
  ## A variable that is endogenous but appears nowhere in the model is a
  ## specification error, not something to ignore.
  expect_error(
    ivsfm(y ~ x1 + x2, ~w1, ~ w1 + w2, d),
    "appear in neither the frontier nor"
  )
  expect_error(
    ivsfm(y ~ x1 + x2, ~x2, ~ x2 + w1, d),
    "both endogenous and an instrument"
  )
  ## Order condition: fewer excluded instruments than endogenous variables.
  expect_error(
    ivsfm(y ~ x1 + x2, ~x2, ~1, d),
    "excluded instrument"
  )
  expect_error(
    ivsfm(y ~ x1 + x2, ~x2, ~ w1 + w2, d, inefdec = NA),
    "must be TRUE or FALSE"
  )
  expect_error(
    ivsfm(y ~ x1 + x2, ~x2, ~ w1 + nosuchvar, d),
    "not found in `data`"
  )
})

test_that("IVLIML and IVCF agree closely but are not the same estimator", {
  d <- iv_data(9, n = 1500)
  a <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "IVLIML"
  )
  b <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "IVCF"
  )
  ## Both are consistent, so they should land close together...
  expect_equal(unname(a$out["x2", "par"]), unname(b$out["x2", "par"]),
    tolerance = 0.05
  )
  ## ...but IVLIML estimates the reduced form jointly, so it has strictly more
  ## parameters. If these ever match exactly, IVLIML has stopped doing that.
  expect_gt(length(a$start_v), length(b$start_v))
})

## Gap A7. On the Hou, Ramalho and Roseta-Palma (2025) Table 2 design,
## "IVLIML" used to walk to the far corner of its own box on a minority of
## draws: rho pinned at -10/sqrt(101) = -0.99503719 with beta2 near 6 against
## a truth of 0.5, returned as an ordinary fit with no error and no warning.
## 27 of 300 draws at rho = 0.8, 6 of 300 at rho = 0.4, 0 of 300 for "IVCF" on
## the same draws. The cause was the starting point, not the search: t starts
## at 0, and span = pmax(10*abs(start_v), 10) then makes that corner exactly
## as far from the start as the true optimum. Seeding from the IVCF solution
## fixes it (0 of 300 in all three designs). These are two of the seeds that
## used to fail.
hou_data <- function(seed, rho, NN = 100L, TT = 20L) {
  set.seed(10000L * as.integer(100 * rho) + seed)
  n <- NN * TT
  x1 <- rnorm(n); z <- rnorm(n); eps <- rnorm(n)
  x2 <- z + eps
  v <- rho * eps + sqrt(1 - rho^2) * rnorm(n)
  u <- rep(abs(rnorm(NN)), each = TT)
  data.frame(y = 0.5 * x1 + 0.5 * x2 + v - u, x1 = x1, x2 = x2, z = z)
}

test_that("IVLIML does not walk into the rho = -10/sqrt(101) corner", {
  rho_bound <- -10 / sqrt(101)
  for (s in c(6L, 25L)) {
    d <- hou_data(s, 0.8)
    f <- ivsfm(y ~ 0 + x1 + x2, endogenous = ~x2, instruments = ~z,
      data = d, model_name = "IVLIML"
    )
    p <- f$out[, "par"]
    ## The signature of the old failure, pinned directly.
    expect_gt(abs(p[["rho_x2"]] - rho_bound), 0.01)
    ## and the consequence of it.
    expect_lt(abs(p[["x2"]] - 0.5), 0.1)
    expect_gt(p[["rho_x2"]], 0)
  }
})
