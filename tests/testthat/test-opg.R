## OPG / BHHH standard errors: keep_objective = TRUE on every entry point whose
## likelihood is separable per observation, and the map from the ESTIMATION
## scale onto the scale coef() reports.
##
## The scale map is the part worth testing hard. ttsfm(), copsfm() and ivsfm()
## estimate LOG sigmas and report sigmas, and ivsfm("IVLIML") estimates
## reduced-form parameters it never reports -- so a score matrix built at
## opt$par is neither on the reported scale nor, in the IVLIML case, of the
## reported length. Before 1.2.1 vcov() ignored that: it named the
## estimation-scale matrix with the reported names and returned it, so
## sqrt(diag(vcov(f))) came back a factor of 1/sigma_u away from the sigma_u
## standard error the same fit printed, and IVLIML failed outright in
## dimnames(). Both are regression-tested below.

cs_data <- function(seed = 1, n = 300) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  data.frame(y = 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n, 0, 0.5) -
    abs(rnorm(n, 0, 1)), x1 = x1, x2 = x2)
}

tt_data <- function(seed = 11, n = 300) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  data.frame(y = 1 + 0.8 * x1 - 0.4 * x2 + rnorm(n, 0, 0.3) -
    abs(rnorm(n, 0, 0.9)) + abs(rnorm(n, 0, 0.7)), x1 = x1, x2 = x2)
}

iv_data <- function(seed = 3, n = 500) {
  set.seed(seed)
  w1 <- rnorm(n); w2 <- rnorm(n); e <- rnorm(n)
  x2 <- 0.7 * w1 - 0.3 * w2 + e
  x1 <- rnorm(n)
  data.frame(y = 0.5 + 0.8 * x1 - 0.6 * x2 + (0.6 * e + rnorm(n, 0, 0.4)) -
    abs(rnorm(n, 0, 1)), x1 = x1, x2 = x2, w1 = w1, w2 = w2)
}

lcm_data <- function(seed = 2, n = 500) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); cl <- rbinom(n, 1, 0.5)
  y <- ifelse(cl == 1, 2 + 1.5 * x1 - 0.5 * x2, -2 + 0.2 * x1 + 1.0 * x2) +
    rnorm(n, 0, 0.4) - abs(rnorm(n, 0, 0.8))
  data.frame(y = y, x1 = x1, x2 = x2)
}

## One fitted example per newly-supported entry point, fitted once and reused:
## several of these are simulated-ML and not cheap.
opg_fits <- local({
  d <- cs_data(); dt <- tt_data(); di <- iv_data()
  suppressWarnings(list(
    zsfm   = zsfm(y ~ x1 + x2, "ZISF", d, keep_objective = TRUE),
    lcsfm  = lcsfm(y ~ x1 + x2, "LCM", lcm_data(), n_class = 2,
      keep_objective = TRUE),
    ttne   = ttsfm(y ~ x1 + x2, "TTNE", dt, keep_objective = TRUE),
    tthn   = ttsfm(y ~ x1 + x2, "TTHN", dt, keep_objective = TRUE),
    copsfm = copsfm(y ~ x1 + x2, data = d, copula = "gaussian",
      keep_objective = TRUE),
    ivcf   = ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = di, model_name = "IVCF", keep_objective = TRUE),
    ivliml = ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = di, model_name = "IVLIML", keep_objective = TRUE)
  ))
})

test_that("keep_objective = TRUE gives every supported entry point a score matrix", {
  for (nm in names(opg_fits)) {
    f <- opg_fits[[nm]]
    G <- estfun.sfareg(f)
    expect_true(is.matrix(G), info = nm)
    expect_true(all(is.finite(G)), info = nm)
    ## One row per observation, one column per ESTIMATED parameter -- which is
    ## not the reported count for IVLIML.
    expect_identical(nrow(G), as.integer(nobs(f)), info = nm)
    expect_identical(ncol(G), length(f$opt$par), info = nm)
  }
})

test_that("the scores satisfy the first-order conditions at the optimum", {
  ## The substantive check that estfun() differences the right thing: at a
  ## maximum the summed score is zero. Scaled by n so the tolerance means the
  ## same thing whatever the sample size.
  for (nm in names(opg_fits)) {
    f <- opg_fits[[nm]]
    G <- estfun.sfareg(f)
    expect_lt(max(abs(colSums(G))) / nrow(G), 1e-3)
  }
})

test_that("vcov(type = 'bhhh') is a usable covariance on the REPORTED scale", {
  for (nm in names(opg_fits)) {
    f <- opg_fits[[nm]]
    V <- vcov(f, type = "bhhh")
    p <- length(coef(f))
    expect_identical(dim(V), c(p, p), info = nm)
    expect_identical(colnames(V), names(coef(f)), info = nm)
    expect_true(all(is.finite(V)), info = nm)
    expect_equal(V, t(V), tolerance = 1e-10)
    expect_true(all(diag(V) > 0), info = nm)
  }
})

test_that("BHHH and Hessian agree on the well-identified slope coefficients", {
  ## The two estimators of the same covariance are free to disagree -- that is
  ## the point of having both, and they disagree most exactly where the
  ## information-matrix equality is weakest. On these fits that is the variance
  ## and mixing parameters: copsfm()'s rho sits on its boundary at -0.95, where
  ## the Hessian standard error is not meaningful and the OPG is 22x larger,
  ## and zsfm()'s gamma implies P(fully efficient) ~ 0.08, so the
  ## zero-inefficiency component is barely identified. Both are real, and
  ## neither is a defect to be tested away.
  ##
  ## The frontier slopes are the parameters that ARE well identified in all of
  ## these designs, and the ones users interpret, so that is where agreement is
  ## genuinely expected.
  for (nm in names(opg_fits)) {
    f <- opg_fits[[nm]]
    ## lcsfm() reports one slope per class ("x1_class1", ...), so match by
    ## prefix rather than by exact name.
    b <- grep("^x[12]", names(coef(f)), value = TRUE)
    r <- (sqrt(diag(vcov(f, type = "bhhh"))) / f$std.errors)[b]
    r <- r[is.finite(r)]
    expect_true(length(r) > 0, info = nm)
    expect_true(all(r > 0.5 & r < 2), info = nm)
  }
})

test_that("the stored Jacobian is the exp() factor it is supposed to be", {
  ## Pins the delta-method factor itself rather than inferring it from a ratio.
  ## These entry points estimate log sigma and report sigma, so
  ## d(reported)/d(estimated) is the reported sigma -- and the Jacobian entry
  ## must equal the coefficient, exactly.
  f <- opg_fits$copsfm
  i <- match(c("sigma_u", "sigma_v"), names(coef(f)))
  expect_equal(unname(f$par_scale[i]), unname(coef(f)[i]), tolerance = 1e-8)
  ## Slopes are reported as estimated, so their factor is exactly 1.
  expect_equal(unname(f$par_scale[match(c("x1", "x2"), names(coef(f)))]),
    c(1, 1), tolerance = 1e-12
  )
  ## ttsfm() under the default "sd" link is the same exp() map on all three.
  g <- opg_fits$tthn
  j <- match(c("sigma_v", "sigma_u", "sigma_w"), names(coef(g)))
  expect_equal(unname(g$par_scale[j]), unname(coef(g)[j]), tolerance = 1e-8)
})

test_that("vcov() agrees with the standard errors the fit reports", {
  ## The regression test for the scale bug. ttsfm(), copsfm() and ivsfm()
  ## estimate log sigmas; vcov() used to return the log-scale matrix under the
  ## reported names, so this ratio came back as 1/sigma for exactly the
  ## transformed parameters and 1 for the rest -- wrong in a way that looks
  ## entirely plausible.
  for (nm in names(opg_fits)) {
    f <- opg_fits[[nm]]
    se_v <- sqrt(diag(vcov(f)))
    keep <- is.finite(se_v) & is.finite(f$std.errors)
    expect_equal(unname(se_v[keep]), unname(f$std.errors[keep]),
      tolerance = 1e-6, info = nm
    )
  }
})

test_that("ivsfm('IVLIML') marginalises its reduced-form parameters", {
  f <- opg_fits$ivliml
  ## 11 estimated, 6 reported: the Jacobian is rectangular and vcov() must
  ## come back square on the reported count, not the estimated one.
  expect_gt(length(f$opt$par), length(coef(f)))
  expect_true(is.matrix(f$par_scale))
  expect_identical(dim(f$par_scale),
    c(length(coef(f)), length(f$opt$par))
  )
  expect_identical(dim(vcov(f)), c(6L, 6L))
  expect_identical(dim(vcov(f, type = "bhhh")), c(6L, 6L))
  ## Marginalising must not be confused with conditioning: taking the block of
  ## the inverse is strictly more uncertain than inverting the block, so the
  ## naive version would understate these.
  G <- estfun.sfareg(f)
  cond <- sqrt(diag(solve(crossprod(G[, seq_len(3L), drop = FALSE]))))
  marg <- sqrt(diag(vcov(f, type = "bhhh")))[seq_len(3L)]
  expect_true(all(marg >= cond * (1 - 1e-8)))
})

test_that("influence_sfa() pairs scores with a covariance on the same scale", {
  for (nm in c("copsfm", "tthn", "ivliml")) {
    inf <- influence_sfa(opg_fits[[nm]])
    expect_identical(ncol(inf$influence), length(coef(opg_fits[[nm]])))
    expect_identical(colnames(inf$influence), names(coef(opg_fits[[nm]])))
  }
})

test_that("without keep_objective the score-based path refuses clearly", {
  f <- suppressWarnings(zsfm(y ~ x1 + x2, "ZISF", cs_data()))
  expect_null(f$objective)
  expect_error(vcov(f, type = "bhhh"), "keep_objective")
  expect_error(estfun.sfareg(f), "keep_objective")
})

test_that("the estimators with no likelihood say so instead of pretending", {
  ## TTNLS is nonlinear least squares and C2SLS is 2SLS with a corrected
  ## intercept; neither maximises a log-likelihood, so neither has a score.
  expect_warning(
    ttsfm(y ~ x1 + x2, "TTNLS", tt_data(), keep_objective = TRUE),
    "no effect"
  )
  expect_warning(
    ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
      data = iv_data(), model_name = "C2SLS", keep_objective = TRUE),
    "no effect"
  )
})
