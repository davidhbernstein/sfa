## Generalized beta of the second kind as the inefficiency distribution:
## sfm(model_name = "NGB2"), gap L20.
##
## Makiela and Mazur (2022), JPA 58:35-54, section 3. The reason for putting
## four parameters on ONE error component is that its limits are most of the
## rest of the package -- half-normal, exponential, gamma, Weibull, half-t --
## so a single fit says which of them the data want.

test_that("the GB2 density integrates to 1 and its quantile inverts its CDF", {
  ld <- sfa:::.gb2_ld; qf <- sfa:::.gb2_q; pf <- sfa:::.gb2_p
  pp <- c(1e-8, 1e-6, 0.05, 0.5, 0.9, 1 - 1e-8, 1 - 1e-12)
  for (nu in c(0.4, 0.8, 3, 50, 1e5)) {
    for (ps in c(0.7, 1, 2)) {
      for (ta in c(0.4, 1, 2.5)) {
        tag <- paste(nu, ps, ta)
        expect_equal(
          integrate(function(u) exp(ld(u, 1.3, nu, ps, ta)), 0, Inf,
            rel.tol = 1e-10)$value,
          1, tolerance = 1e-6, info = tag)
        x <- qf(pp, 1.3, nu, ps, ta)
        ## nu below 1 is a tail so heavy that a naive quantile returns Inf long
        ## before the quantile actually is infinite. A non-finite draw reaches
        ## the likelihood as "steer away", which would make small nu -- the
        ## only region the parameter exists for -- unreachable.
        expect_true(all(is.finite(x)), info = tag)
        expect_equal(pf(x, 1.3, nu, ps, ta), pp, tolerance = 1e-6, info = tag)
      }
    }
  }
  ## The CDF really is the integral of the density, not merely a function that
  ## round-trips with the quantile.
  for (nu in c(1.5, 20)) {
    for (ps in c(0.8, 2)) {
      expect_equal(pf(0.9, 1.3, nu, ps, 0.5),
        integrate(function(u) exp(ld(u, 1.3, nu, ps, 0.5)), 0, 0.9,
          rel.tol = 1e-11)$value, tolerance = 1e-8)
    }
  }
})

test_that("the quantile keeps the dimensions of a matrix of draws", {
  ## REGRESSION. The first version allocated its working vectors with
  ## numeric(length(p)), which drops dim. The simulated-ML path hands this a
  ## matrix; cbind() then recycled the flattened vector against the other
  ## proposal, and the importance weights were formed from misaligned draws --
  ## no error, no warning at the top level, just a wrong likelihood.
  P <- matrix(seq(0.01, 0.99, length.out = 60), nrow = 12)
  q <- sfa:::.gb2_q(P, 1, 30, 2, 1)
  expect_identical(dim(q), dim(P))
  expect_equal(as.numeric(q), sfa:::.gb2_q(as.numeric(P), 1, 30, 2, 1))
})

test_that("GB2 nests the distributions the paper says it nests", {
  ## Their Eqs. (5), (6) and (9): nu -> Inf is the generalized gamma, and the
  ## familiar cases sit inside that. These are LIMITS, not identities, so the
  ## test is that the error falls like 1/nu -- a fixed tolerance would pass on
  ## a density that happened to be close and wrong.
  ld <- sfa:::.gb2_ld
  u <- seq(0.05, 6, by = 0.05); s <- 1.1
  lim <- list(
    "half-normal" = list(c(2, 1), log(2) + dnorm(u, 0, s, log = TRUE), s),
    "exponential" = list(c(1, 1), dexp(u, rate = 1 / s, log = TRUE), s),
    "gamma"       = list(c(1, 3), dgamma(u, shape = 3, scale = s, log = TRUE), s),
    ## GG with tau = psi is Weibull, at scale sigma psi^(1/psi).
    "Weibull"     = list(c(1.7, 1.7), dweibull(u, shape = 1.7, scale = s, log = TRUE),
                         s * 1.7^(-1 / 1.7))
  )
  for (nm in names(lim)) {
    sh <- lim[[nm]][[1]]; want <- lim[[nm]][[2]]; sg <- lim[[nm]][[3]]
    e4 <- max(abs(ld(u, sg, 1e4, sh[1], sh[2]) - want))
    e6 <- max(abs(ld(u, sg, 1e6, sh[1], sh[2]) - want))
    expect_lt(e6, 1e-3, label = nm)
    ## First order in 1/nu: a hundredfold nu buys a hundredfold accuracy.
    expect_equal(e4 / e6, 100, tolerance = 0.05, info = nm)
  }
  ## The half-Student t is EXACT, not a limit: tau = 1, psi = 2, nu = df.
  expect_equal(ld(u, s, 4, 2, 1), log(2) + dt(u / s, df = 4, log = TRUE) - log(s),
    tolerance = 1e-12)
})

test_that("sfm(model_name = 'NGB2') returns a well-formed fit", {
  skip_on_cran()
  set.seed(5)
  n <- 400
  x1 <- rnorm(n); x2 <- rnorm(n)
  d <- data.frame(y = 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n, 0, 0.4) -
      abs(rnorm(n, 0, 1)), x1 = x1, x2 = x2)
  f <- suppressWarnings(sfm(y ~ x1 + x2, data = d, model_name = "NGB2",
    maxit.bobyqa = 400, maxit.optim = 200))

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "NGB2")
  expect_identical(rownames(f$out),
    c("sigv", "sigu", "nu", "psi", "tau", "(Intercept)", "x1", "x2"))
  expect_true(all(f$out[c("sigv", "sigu", "nu", "psi", "tau"), "par"] > 0))
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_equal(length(f$exp_u_hat), n)
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_equal(nobs(f), n)
  ## Five non-frontier parameters precede the slopes; if the offset were wrong
  ## the frontier would be fitted with shape parameters in place of betas.
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.15)
  expect_equal(unname(f$out["x2", "par"]), -0.4, tolerance = 0.15)
})

test_that("a large fitted nu is reported as the generalized-gamma limit", {
  ## nu has no interior maximum once the tail stops being heavy: past a few
  ## hundred the density has stopped moving, so the number is a statement
  ## about the family and not an estimate. Warned about, as tHN's nu is.
  skip_on_cran()
  set.seed(6)
  n <- 300
  x1 <- rnorm(n)
  d <- data.frame(y = 0.5 + 0.8 * x1 + rnorm(n, 0, 0.4) - abs(rnorm(n, 0, 1)),
    x1 = x1)
  w <- tryCatch({
    f <- sfm(y ~ x1, data = d, model_name = "NGB2",
      maxit.bobyqa = 300, maxit.optim = 150)
    NA_character_
  }, warning = function(x) conditionMessage(x))
  ## Either it warned about nu, or nu came back small enough not to warrant it.
  if (!is.na(w) && grepl("nu converged", w)) {
    expect_match(w, "generalized-gamma limit")
  } else {
    f <- suppressWarnings(sfm(y ~ x1, data = d, model_name = "NGB2",
      maxit.bobyqa = 300, maxit.optim = 150))
    expect_lte(unname(f$out["nu", "par"]), 200)
  }
})
