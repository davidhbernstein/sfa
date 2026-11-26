## Discrete determinants get discrete changes, not derivatives. A dummy's
## effect is E[u | z_k = 1] - E[u | z_k = 0] per observation, other covariates
## at their observed values; a named ordinal gets the change per unit across
## its range; a factor's dummies are switched as a group against the
## reference. The fits below are hand-built: marginal_effects() reads only the
## stored design blocks, so the arithmetic is pinned against a formula written
## out here rather than against wherever an optimizer lands.

set.seed(11)
n_me <- 60
z_c <- rnorm(n_me)
z_d <- rbinom(n_me, 1, 0.4)
z_o <- sample(1:3, n_me, TRUE)
Zm <- cbind("(Intercept)" = 1, zc = z_c, zd = z_d, zo = z_o)

fake_fit <- function(z_spec, ...) {
  structure(c(list(z_spec = z_spec, model_name = "fake"), list(...)), class = "sfareg")
}
spec <- function(Z, delta, link = "sd", family = "halfnormal") {
  list(Z = Z, delta = stats::setNames(delta, colnames(Z)), link = link, family = family)
}
set_col <- function(Z, k, v) {
  Z[, k] <- v
  Z
}
## E[u] and Var[u] of the half-normal and exponential at a design, written out.
mom_hn_e <- function(Z, d, link) {
  s <- if (link == "sd") exp(drop(Z %*% d)) else sqrt(exp(drop(Z %*% d)))
  list(hn = list(e = s * sqrt(2 / pi), v = s^2 * (1 - 2 / pi)), ex = list(e = s, v = s^2))
}

test_that("a dummy is a first difference, for both families and both links", {
  d <- c(0.2, 0.3, -0.5, 0.1)
  for (link in c("sd", "var")) for (fam in c("halfnormal", "exponential")) {
    me <- marginal_effects(fake_fit(spec(Zm, d, link, fam)))
    key <- if (fam == "halfnormal") "hn" else "ex"
    hi <- mom_hn_e(set_col(Zm, "zd", 1), d, link)[[key]]
    lo <- mom_hn_e(set_col(Zm, "zd", 0), d, link)[[key]]
    expect_equal(me[["DeltaE_u.zd"]], hi$e - lo$e, tolerance = 1e-13, info = paste(link, fam))
    expect_equal(me[["DeltaVar_u.zd"]], hi$v - lo$v, tolerance = 1e-13, info = paste(link, fam))
    expect_equal(unname(attr(me, "average")["DeltaE_u.zd"]), mean(hi$e - lo$e),
                 tolerance = 1e-13)
    expect_identical(unname(attr(me, "effect_type")["DeltaE_u.zd"]), "discrete")
  }
})

test_that("continuous columns are exactly what they were, and discrete = NULL is the old output", {
  d <- c(0.2, 0.3, -0.5, 0.1)
  f <- fake_fit(spec(Zm, d))
  new <- marginal_effects(f)
  old <- marginal_effects(f, discrete = NULL)
  expect_identical(new[["dE_u.dzc"]], old[["dE_u.dzc"]])
  expect_identical(new[["dVar_u.dzo"]], old[["dVar_u.dzo"]])
  ## The old output differentiates the dummy, with the derivative's name.
  expect_true("dE_u.dzd" %in% names(old))
  expect_false(any(grepl("^Delta", names(old))))
  expect_true(all(attr(old, "effect_type") == "derivative"))
  expect_identical(marginal_effects(f, discrete = FALSE), old)
})

test_that("a named ordinal is the change per unit across its range", {
  d <- c(0.2, 0.3, -0.5, 0.1)
  for (link in c("sd", "var")) {
    me <- marginal_effects(fake_fit(spec(Zm, d, link)), discrete = c("auto", "zo"))
    hi <- mom_hn_e(set_col(Zm, "zo", 3), d, link)$hn
    lo <- mom_hn_e(set_col(Zm, "zo", 1), d, link)$hn
    ## The erratum's definition: mean of [E(Age = 3) - E(Age = 1)] / 2.
    expect_equal(unname(attr(me, "average")["DeltaE_u.zo"]), mean((hi$e - lo$e) / 2),
                 tolerance = 1e-13)
    expect_identical(unname(attr(me, "effect_type")["DeltaE_u.zo"]), "one step")
    ## Naming a column without "auto" makes the 0/1 column a derivative again.
    only <- marginal_effects(fake_fit(spec(Zm, d, link)), discrete = "zo")
    expect_true("dE_u.dzd" %in% names(only))
  }
})

test_that("an unknown name in `discrete` is an error that lists the choices", {
  f <- fake_fit(spec(Zm, c(0.2, 0.3, -0.5, 0.1)))
  expect_error(marginal_effects(f, discrete = "nope"), "not among")
  expect_error(marginal_effects(f, discrete = 1), "must be")
})

## Truncated normal: E[u] = mu + s lambda(mu/s), Var[u] = s^2 (1 - lambda (lambda + a)).
tn_mom <- function(mu, s) {
  a <- mu / s
  lam <- exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
  list(e = mu + s * lam, v = s^2 * (1 - lam * (lam + a)))
}

test_that("truncated normal: mu and sigma_u move together, incl. a dummy only in mu", {
  Zu <- Zm[, c("(Intercept)", "zc", "zd")]
  Zmu <- cbind("(Intercept)" = 1, zd = z_d, zo = z_o)
  du <- c(0.1, 0.2, -0.3)
  dmu <- c(0.4, 0.5, -0.2)
  f <- fake_fit(spec(Zu, du, "sd", "truncnormal"),
                mu_spec = list(Z = Zmu, delta = stats::setNames(dmu, colnames(Zmu))))
  me <- marginal_effects(f)
  at <- function(v) {
    tn_mom(drop(set_col(Zmu, "zd", v) %*% dmu), exp(drop(set_col(Zu, "zd", v) %*% du)))
  }
  expect_equal(me[["DeltaE_u.zd"]], at(1)$e - at(0)$e, tolerance = 1e-12)
  expect_equal(me[["DeltaVar_u.zd"]], at(1)$v - at(0)$v, tolerance = 1e-12)
  ## zo sits in the mu design alone and is still reported.
  me2 <- marginal_effects(f, discrete = c("auto", "zo"))
  m3 <- tn_mom(drop(set_col(Zmu, "zo", 3) %*% dmu), exp(drop(Zu %*% du)))
  m1 <- tn_mom(drop(set_col(Zmu, "zo", 1) %*% dmu), exp(drop(Zu %*% du)))
  expect_equal(me2[["DeltaE_u.zo"]], (m3$e - m1$e) / 2, tolerance = 1e-12)
  ## A dummy in the mu design ONLY (zd removed from the scale).
  f2 <- fake_fit(spec(Zu[, 1:2], du[1:2], "sd", "truncnormal"),
                 mu_spec = list(Z = Zmu, delta = stats::setNames(dmu, colnames(Zmu))))
  me3 <- marginal_effects(f2)
  s0 <- exp(drop(Zu[, 1:2] %*% du[1:2]))
  hi <- tn_mom(drop(set_col(Zmu, "zd", 1) %*% dmu), s0)
  lo <- tn_mom(drop(set_col(Zmu, "zd", 0) %*% dmu), s0)
  expect_equal(me3[["DeltaE_u.zd"]], hi$e - lo$e, tolerance = 1e-12)
})

test_that("scaling: a dummy in the scaling design is a first difference of h E[u*]", {
  Zs <- Zm[, c("zc", "zd"), drop = FALSE]
  ds <- c(0.3, -0.4)
  f <- fake_fit(spec(Zm[, 1, drop = FALSE], 0.1, "sd", "truncnormal"),
                s_spec = list(Z = Zs, delta = stats::setNames(ds, colnames(Zs))),
                mu_spec = list(Z = Zm[, 1, drop = FALSE], delta = c("(Intercept)" = 0.5)))
  me <- marginal_effects(f)
  base <- tn_mom(0.5, exp(0.1))
  h <- function(v) exp(drop(set_col(Zs, "zd", v) %*% ds))
  expect_equal(me[["DeltaE_u.zd"]], (h(1) - h(0)) * base$e, tolerance = 1e-12)
  expect_equal(me[["DeltaVar_u.zd"]], (h(1)^2 - h(0)^2) * base$v, tolerance = 1e-12)
  ## The continuous scaling column is still delta_s * E[u].
  expect_equal(me[["dE_u.dzc"]], 0.3 * me[["E_u"]], tolerance = 1e-12)
})

test_that("component = 'h' gets discrete changes on the persistent block", {
  d <- c(0.2, 0.3, -0.5, 0.1)
  f <- structure(list(z_spec = spec(Zm, d, "var"), z_spec_h = spec(Zm, -d, "var"),
                      model_name = "GTRE_Z"), class = "sfareg")
  mh <- marginal_effects(f, component = "h")
  hi <- mom_hn_e(set_col(Zm, "zd", 1), -d, "var")$hn
  lo <- mom_hn_e(set_col(Zm, "zd", 0), -d, "var")$hn
  expect_equal(mh[["DeltaE_h.zd"]], hi$e - lo$e, tolerance = 1e-13)
  expect_true("DeltaVar_h.zd" %in% names(mh))
})

test_that("a factor's dummies are switched as a group against the reference", {
  fct <- factor(sample(c("a", "b", "c"), n_me, TRUE))
  Zf <- model.matrix(~ z_c + fct)
  colnames(Zf)[2] <- "zc"
  d <- c(0.2, 0.3, -0.4, 0.25)
  me <- marginal_effects(fake_fit(spec(Zf, d)))
  lvl <- function(b, c) {
    Z <- Zf
    Z[, "fctb"] <- b
    Z[, "fctc"] <- c
    unname(mom_hn_e(Z, d, "sd")$hn$e)
  }
  ## Level b against the reference: fctc held at 0, never both switched on.
  expect_equal(me[["DeltaE_u.fctb"]], lvl(1, 0) - lvl(0, 0), tolerance = 1e-13)
  expect_equal(me[["DeltaE_u.fctc"]], lvl(0, 1) - lvl(0, 0), tolerance = 1e-13)
  expect_identical(unname(attr(me, "effect_type")["DeltaE_u.fctb"]), "factor level")
  ## Without the factor attributes (a pipe segment's design) each dummy is
  ## switched alone, as the erratum's calculation did.
  Zp <- Zf
  attr(Zp, "assign") <- NULL
  attr(Zp, "contrasts") <- NULL
  mp <- marginal_effects(fake_fit(spec(Zp, d)))
  single <- mom_hn_e(set_col(Zp, "fctb", 1), d, "sd")$hn$e -
    mom_hn_e(set_col(Zp, "fctb", 0), d, "sd")$hn$e
  expect_equal(mp[["DeltaE_u.fctb"]], unname(single), tolerance = 1e-13)
  expect_identical(unname(attr(mp, "effect_type")["DeltaE_u.fctb"]), "discrete")
})

test_that("on a real NHN_Z fit a dummy determinant reports a first difference", {
  skip_on_cran()
  dat <- data_gen_cs(N = 300, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  set.seed(1)
  dat$dmy <- rbinom(300, 1, 0.5)
  fit <- suppressWarnings(sfm(y_pcs_z ~ x1 + x2 | z + dmy, data = dat, model_name = "NHN_Z"))
  me <- marginal_effects(fit)
  Z <- fit$z_spec$Z
  dl <- fit$z_spec$delta
  Eu <- function(Z) unname(sqrt(2 / pi) * exp(drop(Z %*% dl)))
  expect_equal(me[["DeltaE_u.dmy"]], Eu(set_col(Z, "dmy", 1)) - Eu(set_col(Z, "dmy", 0)),
               tolerance = 1e-12)
  expect_identical(me[["dE_u.dz"]], marginal_effects(fit, discrete = NULL)[["dE_u.dz"]])
})
