## vcov(type = ) offers four covariances: "hessian", "bhhh", "sandwich" and
## "clustered". The last two are computed in the package rather than left to
## sandwich::sandwich(), because that function composes bread and meat itself
## and knows nothing about par_index/par_scale -- the map from the ESTIMATION
## scale the scores live on to the REPORTED scale coef() returns.
##
## Composing a reported-scale bread with an estimation-scale meat gave finite,
## plausible, WRONG numbers on every model whose two scales differ: on
## ttsfm("TTNE") the sigma standard errors were 29% low, 18% low and 26% HIGH,
## each in a different direction, while the betas (Jacobian 1) were correct.
## On ivsfm("IVLIML"), which estimates 11 and reports 6, it could not conform
## at all and errored. Both are pinned below.

.vt_cs <- function(n = 300L, seed = 5L) {
  as.data.frame(data_gen_cs(N = n, rand = seed, sig_u = 1, sig_v = 0.5,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1))
}

test_that("on an identity-scale fit the new types match the sandwich package exactly", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = .vt_cs(), keep_objective = TRUE))
  g <- rep_len(seq_len(10L), nobs(f))

  ## sfm() reports what it estimates, so package and sandwich must agree to
  ## machine precision. This is the anchor: an independent implementation of
  ## the same quantity.
  expect_equal(unname(vcov(f, type = "sandwich")),
    unname(sandwich::sandwich(f)), tolerance = 1e-10)
  ## The cluster adjustment is M/(M-1) only -- vcovCL() does NOT also apply the
  ## HC1 factor (n-1)/(n-k). Getting that wrong left the two disagreeing by
  ## exactly (n-1)/(n-k), so this assertion is what pins the correction.
  expect_equal(unname(vcov(f, type = "clustered", cluster = g)),
    unname(sandwich::vcovCL(f, cluster = g)), tolerance = 1e-10)
})

test_that("on a fit whose scales differ, sandwich is delta-mapped to the reported scale", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  f <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2, model_name = "TTNE",
    data = .vt_cs(), keep_objective = TRUE))

  G <- sandwich::estfun(f)
  J <- as.numeric(attr(G, "par_scale"))
  expect_true(any(abs(J - 1) > 1e-8))          # otherwise this tests nothing
  Hi <- solve(f$opt$hessian)
  hand <- diag(J) %*% (Hi %*% crossprod(G) %*% Hi) %*% diag(J)

  expect_equal(unname(vcov(f, type = "sandwich")), unname(hand), tolerance = 1e-8)

  ## Teeth: the mixed-scale product the old bread() produced is a DIFFERENT
  ## matrix, so this fails loudly if the mapping is ever dropped.
  mixed <- (vcov(f) * nobs(f)) %*% (crossprod(G) / nobs(f)) %*% (vcov(f) * nobs(f)) / nobs(f)
  expect_false(isTRUE(all.equal(unname(vcov(f, type = "sandwich")), unname(mixed),
    tolerance = 1e-6)))
})

test_that("ivsfm('IVLIML') reports 6 parameters from 11 estimated, on every type", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  ## Same fixture the endogeneity tests use, so the fit is one the suite
  ## already exercises rather than one invented here.
  set.seed(9)
  n <- 800L
  z <- stats::rnorm(n); x1 <- stats::rnorm(n); eps <- stats::rnorm(n)
  x2 <- z + eps
  v <- 0.6 * eps + sqrt(1 - 0.6^2) * stats::rnorm(n)
  u <- abs(stats::rnorm(n))
  d <- data.frame(y = 0.5 * x1 + 0.5 * x2 + v - u, x1 = x1, x2 = x2, z = z)
  f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~z,
    model_name = "IVLIML", data = d, keep_objective = TRUE))
  skip_if(is.null(f$opt) || is.null(f$objective))

  p <- length(coef(f))
  G <- sandwich::estfun(f)
  skip_if(ncol(G) == p)                        # only meaningful when they differ

  for (ty in c("hessian", "bhhh", "sandwich")) {
    V <- vcov(f, type = ty)
    expect_equal(dim(V), c(p, p), info = ty)
    expect_equal(colnames(V), names(coef(f)), info = ty)
  }
  V <- vcov(f, type = "clustered", cluster = rep_len(seq_len(6L), nrow(G)))
  expect_equal(dim(V), c(p, p))

  ## sandwich::sandwich() is now internally consistent -- estimation scale on
  ## both sides -- where before the dimensions did not conform at all.
  expect_equal(dim(sandwich::sandwich(f)), c(ncol(G), ncol(G)))
})

test_that("every type returns a symmetric p x p matrix named as coef()", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = .vt_cs(), keep_objective = TRUE))
  p <- length(coef(f))
  g <- rep_len(seq_len(8L), nobs(f))
  for (ty in c("hessian", "bhhh", "sandwich", "clustered")) {
    V <- if (identical(ty, "clustered")) vcov(f, type = ty, cluster = g) else vcov(f, type = ty)
    expect_equal(dim(V), c(p, p), info = ty)
    expect_equal(rownames(V), names(coef(f)), info = ty)
    expect_equal(colnames(V), names(coef(f)), info = ty)
    expect_equal(V, t(V), tolerance = 1e-10, info = ty)
    expect_true(all(is.finite(sqrt(diag(V)))), info = ty)
  }
})

test_that("clustered refuses the mistakes that are easy to make", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, data = .vt_cs(), keep_objective = TRUE))
  expect_error(vcov(f, type = "clustered"), "`cluster` is required")
  expect_error(vcov(f, type = "clustered", cluster = seq_len(7L)), "but the score matrix has")
  expect_error(vcov(f, type = "clustered", cluster = rep(1L, nobs(f))), "single group")
})

test_that("a panel fit clusters over FIRMS, and says so when handed N*T", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  d <- as.data.frame(data_gen_p(t = 5, N = 40, rand = 5, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  f <- suppressWarnings(psfm(y_tre ~ x1 + x2, model_name = "TRE", data = d,
    individual = "name", keep_objective = TRUE))
  n_firms <- nrow(sandwich::estfun(f))
  expect_equal(n_firms, 40L)

  ## The natural thing to reach for is one row per firm-YEAR, and it is the
  ## wrong length by construction. The message has to name the unit.
  expect_error(vcov(f, type = "clustered", cluster = rep_len(seq_len(4L), nrow(d))),
    "FIRM")
  V <- vcov(f, type = "clustered", cluster = rep_len(seq_len(4L), n_firms))
  expect_equal(dim(V), c(length(coef(f)), length(coef(f))))
})
