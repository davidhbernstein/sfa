## Copula stochastic frontier: copsfm(), gap K5.

cop_gen <- function(seed, n, rho, su = 1, sv = 0.4, b = c(0.5, 0.8, -0.4)) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ## Gaussian copula between the marginals: draw correlated normals, push each
  ## through its own marginal quantile function. w1 = Phi(z1) is F_V(v) and
  ## w2 = Phi(z2) is F_U(u), so (w1, w2) carry exactly the dependence rho.
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  v <- sv * z1
  u <- su * qnorm((1 + pnorm(z2)) / 2)
  data.frame(y = b[1] + b[2] * x1 + b[3] * x2 + v - u, x1 = x1, x2 = x2)
}

test_that("each copula density integrates to 1 and is 1 at independence", {
  ## The check that catches a mistranscribed density. Both are asserted for
  ## every family, which is why families whose density could not be verified
  ## are deliberately absent from .cop_logc().
  gl <- .gauss_legendre_01(80L)
  G <- expand.grid(a = gl$nodes, b = gl$nodes)
  W <- as.vector(outer(gl$weights, gl$weights))
  for (fam in c("gaussian", "fgm")) {
    for (p in c(-0.6, -0.2, 0, 0.3, 0.7)) {
      I <- sum(W * exp(.cop_logc(G$a, G$b, p, fam)))
      expect_equal(I, 1, tolerance = 1e-4, info = paste(fam, p))
    }
    ## At the independence parameter the density is exactly 1, so log c is 0.
    expect_equal(max(abs(.cop_logc(G$a, G$b, 0, fam))), 0, info = fam)
  }
})

test_that("the composed density reduces to normal/half-normal at independence", {
  ## With c = 1 the quadrature must reproduce the closed form the rest of the
  ## package uses. This checks the composition and the change of variable
  ## u = t/(1-t) together.
  gl <- .gauss_legendre_01(128L)
  tt <- gl$nodes; wt <- gl$weights
  uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  su <- 1; sv <- 0.4; e <- seq(-4, 2, length.out = 9)

  U <- matrix(uu, length(e), length(uu), byrow = TRUE)
  V <- e + U
  lg <- dnorm(V, 0, sv, log = TRUE) + log(2) + dnorm(U, 0, su, log = TRUE)
  lg <- sweep(lg, 2, log(wt) + log(jac), "+")
  got <- .log_row_sum_exp(lg)

  sig <- sqrt(su^2 + sv^2); lam <- su / sv
  want <- log(2) - log(sig) + dnorm(e / sig, log = TRUE) +
    pnorm(-e * lam / sig, log.p = TRUE)
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("copsfm returns a well-formed sfareg", {
  skip_on_cran()
  d <- cop_gen(4, 600, 0.5)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "COP")
  expect_identical(f$copula, "gaussian")
  expect_equal(ncol(f$out), 3L)
  expect_identical(rownames(f$out),
    c("(Intercept)", "x1", "x2", "sigma_u", "sigma_v", "rho"))
  expect_gt(f$out["sigma_u", "par"], 0)
  expect_lt(abs(f$out["rho", "par"]), 1)
  expect_equal(length(f$jlms), nrow(d))
  expect_equal(f$efficiency, exp(-f$jlms))
  expect_equal(nobs(f), nrow(d))
  ## fgm names its parameter differently, and reports it.
  g <- copsfm(y ~ x1 + x2, data = d, copula = "fgm")
  expect_true("theta" %in% rownames(g$out))
})

test_that("the frontier slopes are recovered even where rho is not", {
  skip_on_cran()
  ## The dependence parameter needs a large sample (see ?copsfm); the SLOPES do
  ## not, and separating the two is the point of this test.
  d <- cop_gen(4, 2000, 0.6)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.08)
  expect_equal(unname(f$out["x2", "par"]), -0.4, tolerance = 0.08)
})

test_that("copsfm rejects malformed calls", {
  d <- cop_gen(2, 300, 0.3)
  expect_error(copsfm(~ x1 + x2, data = d), "two-sided formula")
  expect_error(copsfm(y ~ x1 | x2, data = d), "must not contain a `\\|` segment")
  expect_error(copsfm(y ~ x1 + x2, data = as.matrix(d)), "must be a data.frame")
  expect_error(copsfm(y ~ x1 + x2, data = d, inefdec = NA), "must be TRUE or FALSE")
  expect_error(copsfm(y ~ x1 + x2, data = d, n_nodes = 4), ">= 16")
  ## "clayton" used to belong here, as a family deliberately not implemented.
  ## It is implemented now, so the guard moves to one that still is not: the
  ## point is that an unrecognised family is refused, not that this particular
  ## one is missing.
  expect_error(copsfm(y ~ x1 + x2, data = d, copula = "bb1"))
  expect_error(copsfm(y ~ x1 + x2, data = d, copula = "tawn"))
})

## ---------------------------------------------------------------------------
## Frank, Clayton, Gumbel, Joe and the rotations (added 2026-09-04)
## ---------------------------------------------------------------------------
##
## The families the header of .cop_logc() used to list as deliberately absent.
## They are admitted now because each is checked against something OTHER than
## itself: its own CDF. A density that integrates to 1 can still be the wrong
## density; a density that equals the second mixed partial of the right CDF at
## forty points cannot be.

test_that("each new density IS the second mixed partial of its own CDF", {
  CDF <- list(
    frank = function(u, v, th)
      -1 / th * log1p(expm1(-th * u) * expm1(-th * v) / expm1(-th)),
    clayton = function(u, v, th) (u^(-th) + v^(-th) - 1)^(-1 / th),
    gumbel = function(u, v, th) exp(-(((-log(u))^th + (-log(v))^th)^(1 / th))),
    joe = function(u, v, th)
      1 - ((1 - u)^th + (1 - v)^th - (1 - u)^th * (1 - v)^th)^(1 / th)
  )
  num <- function(f, u, v, th, h = 1e-5) {
    (f(u + h, v + h, th) - f(u + h, v - h, th) -
       f(u - h, v + h, th) + f(u - h, v - h, th)) / (4 * h^2)
  }
  pars <- list(frank = c(-8, -2, 2, 8), clayton = c(0.3, 1, 3, 8),
    gumbel = c(1.3, 2, 4), joe = c(1.3, 2, 5))
  g <- expand.grid(u = seq(0.1, 0.9, by = 0.2), v = seq(0.1, 0.9, by = 0.2))
  for (fam in names(CDF)) {
    for (th in pars[[fam]]) {
      a <- exp(.cop_logc(g$u, g$v, th, fam))
      b <- mapply(function(u, v) num(CDF[[fam]], u, v, th), g$u, g$v)
      expect_equal(a, b, tolerance = 1e-4, info = paste(fam, th))
    }
  }
})

test_that("the new densities integrate to 1 and are 1 at independence", {
  gl <- .gauss_legendre_01(120L)
  G <- expand.grid(a = gl$nodes, b = gl$nodes)
  W <- as.vector(outer(gl$weights, gl$weights))
  pars <- list(frank = c(-8, -2, 2, 8), clayton = c(0.3, 1), gumbel = c(1.3, 2),
    joe = c(1.3, 2))
  for (fam in names(pars)) {
    for (p in pars[[fam]]) {
      I <- sum(W * exp(.cop_logc(G$a, G$b, p, fam)))
      ## Looser than the Gaussian/FGM tolerance ON PURPOSE: these densities
      ## concentrate in a corner, and what is left is quadrature error rather
      ## than density error -- which is why the mixed-partial test above, and
      ## not this one, is the check that would catch a wrong formula.
      expect_equal(I, 1, tolerance = 5e-3, info = paste(fam, p))
    }
  }
  ind <- c(frank = 0, clayton = 0, gumbel = 1, joe = 1)
  for (fam in names(ind)) {
    expect_equal(max(abs(.cop_logc(G$a, G$b, ind[[fam]], fam))), 0, info = fam)
  }
})

test_that("rotations reverse the sign of dependence, which is why they exist", {
  ## Clayton, Gumbel and Joe carry only POSITIVE dependence. Nothing rules out
  ## a negative association between noise and inefficiency, so without the
  ## rotations those families could only ever report the independence boundary
  ## against negatively dependent data.
  set.seed(1)
  n <- 1e5
  u <- runif(n); v <- runif(n)
  rho_s <- function(fam, th) {
    w <- exp(.cop_logc_rot(u, v, th, fam))
    w[!is.finite(w)] <- 0
    12 * sum(w * (u - 0.5) * (v - 0.5)) / sum(w)
  }
  for (base in c("clayton", "gumbel", "joe")) {
    th <- if (base == "clayton") 3 else 2.5
    pos <- rho_s(base, th)
    expect_gt(pos, 0.3)
    expect_lt(rho_s(paste0(base, "90"), th), -0.3)
    expect_lt(rho_s(paste0(base, "270"), th), -0.3)
    ## The survival rotation preserves the sign rather than flipping it.
    expect_gt(rho_s(paste0(base, "180"), th), 0.3)
  }
})

test_that("every advertised family actually fits", {
  skip_on_cran()
  d <- cop_gen(seed = 4, n = 250, rho = 0.3)
  fams <- eval(formals(copsfm)$copula)
  for (fam in fams) {
    f <- suppressWarnings(copsfm(y ~ x1 + x2, data = d, copula = fam,
      n_nodes = 32))
    expect_s3_class(f, "sfareg")
    expect_true(is.finite(as.numeric(logLik(f))), info = fam)
    if (identical(fam, "independent")) {
      ## No dependence parameter exists to be inside a bound.
      expect_true(is.na(f$copula_par))
      expect_false("rho" %in% rownames(f$out))
    } else {
      sp <- sfa:::.cop_spec(fam)
      expect_gte(f$copula_par, sp$lo)
      expect_lte(f$copula_par, sp$hi)
    }
  }
})

test_that("families that do not recover their parameter warn, and the others do not", {
  skip_on_cran()
  d <- cop_gen(seed = 5, n = 150, rho = 0.3)
  fit <- function(fam) copsfm(y ~ x1 + x2, data = d, copula = fam, n_nodes = 32)

  ## Verified to recover: silent. Asserted with expect_silent-style intent but
  ## via the message text, because copsfm() can emit unrelated optimiser
  ## warnings at small n_nodes and those are not what this test is about.
  for (fam in c("gaussian", "fgm", "frank")) {
    w <- tryCatch({ suppressMessages(fit(fam)); NA_character_ },
      warning = function(x) conditionMessage(x))
    expect_false(isTRUE(grepl("did not recover|has not been measured", w)), info = fam)
  }

  ## Measured to collapse: the warning names the measured rate, so the user is
  ## told how bad it is rather than merely that something is wrong.
  expect_warning(fit("gumbel"), "did not recover.*36%")
  expect_warning(fit("joe"), "did not recover.*40%")
  expect_warning(fit("clayton270"), "did not recover.*56%")
  expect_warning(fit("gumbel90"), "did not recover.*60%")

  ## Not measured: says so, rather than implying either result.
  expect_warning(fit("joe180"), "has not been measured")
  expect_warning(fit("clayton90"), "has not been measured")

  ## A missing name must give NA, not an error -- `[[` on a named vector throws.
  expect_true(is.na(unname(sfa:::.COP_COLLAPSE["gumbel270"])))
  expect_silent(sfa:::.cop_warn_family("frank"))
})

## ---------------------------------------------------------------------------
## Skewed noise and non-half-normal inefficiency (added 2026-09-10, gap L18)
## ---------------------------------------------------------------------------
##
## Bonanno and Domma (2022), Mathematics 10:3876. Their contribution beyond the
## copula is a GENERALIZED LOGISTIC noise, whose own skewness parameter gives a
## positive residual skew somewhere to go other than sigma_u = 0, paired with
## an exponential inefficiency.

test_that("the generalized logistic has the moments the paper claims", {
  ## E[V] = 0 for EVERY (alpha, delta) is the property the whole specification
  ## rests on: without it the noise has a mean and is not separable from the
  ## intercept. The other two identities pin the parameterization -- a shifted
  ## or rescaled GL would still integrate to 1 and would still be unimodal.
  for (a in c(0.4, 1, 2.5)) {
    for (d in c(0.3, 1, 2)) {
      f <- function(v) exp(sfa:::.gl_ld(v, a, d))
      expect_equal(integrate(f, -Inf, Inf)$value, 1, tolerance = 1e-7,
        info = paste(a, d))
      expect_equal(integrate(function(v) v * f(v), -Inf, Inf)$value, 0,
        tolerance = 1e-7, info = paste(a, d))
      expect_equal(integrate(function(v) v^2 * f(v), -Inf, Inf)$value,
        d^2 * (trigamma(a) + trigamma(1)), tolerance = 1e-6, info = paste(a, d))
      expect_equal(integrate(function(v) v^3 * f(v), -Inf, Inf)$value,
        d^3 * (psigamma(a, 2) - psigamma(1, 2)), tolerance = 1e-6,
        info = paste(a, d))
    }
  }
  ## alpha = 1 is the ordinary logistic, which is what vdist = "logistic" is.
  v <- seq(-6, 6, by = 0.5)
  expect_equal(sfa:::.gl_ld(v, 1, 0.8), dlogis(v, 0, 0.8, log = TRUE))
  ## The quantile inverts the CDF, including in both tails.
  p <- c(1e-8, 0.01, 0.3, 0.5, 0.9, 1 - 1e-8)
  expect_equal(exp(sfa:::.gl_lp(sfa:::.gl_q(p, 2.5, 0.7), 2.5, 0.7)), p,
    tolerance = 1e-10)
  ## The sign of the skew follows alpha, which is the point of the parameter.
  expect_lt(psigamma(0.5, 2) - psigamma(1, 2), 0)
  expect_gt(psigamma(3, 2) - psigamma(1, 2), 0)
})

## The package's composed density, assembled from the same pieces copsfm()'s
## internal closure uses. Kept here rather than reaching into the closure so
## that the check is against the PARTS the estimator is built from.
.cop_test_dens <- function(e, av, dv, du, th, S, nn = 256) {
  g <- sfa:::.gauss_legendre_01(as.integer(nn))
  tt <- g$nodes; wt <- g$weights; uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  U <- matrix(uu, length(e), length(uu), byrow = TRUE)
  eps <- S * e                          # e is the observed y - Xb
  vv <- S * (eps + U)
  lg <- sfa:::.cop_v_ld(vv, dv, av, "glogistic") +
    sfa:::.cop_u_ld(U, du, "exponential")
  w1 <- sfa:::.cop_v_p(vv, dv, av, "glogistic")
  w2 <- sfa:::.cop_u_p(U, du, "exponential")
  lg <- lg + matrix(sfa:::.cop_logc_rot(as.numeric(w1), as.numeric(w2), th, "fgm"),
    nrow = length(e))
  exp(sfa:::.log_row_sum_exp(sweep(lg, 2, log(wt) + log(jac), "+")))
}

test_that("the composed density IS Bonanno and Domma's Theorems 1 and 2", {
  ## The strongest check available for this model: the paper publishes the
  ## density in closed form, as a sum of four Gauss hypergeometric terms, and
  ## the package computes it by quadrature. Neither derives from the other.
  ## gsl::hyperg_2F1 needs |x| < 1 and here the argument is -k1, which is
  ## unbounded, so Pfaff's transformation moves it into (0, 1).
  F21 <- function(a, b, cc, s) (1 - s)^(-b) * gsl::hyperg_2F1(cc - a, b, cc, s / (s - 1))
  k1f <- function(e, av, dv) exp(-(e + dv * (digamma(av) - digamma(1))) / dv)

  bd_prod <- function(e, av, dv, du, th) {              # Theorem 1, eps = v - u
    k <- k1f(e, av, dv); r <- dv / du
    T <- function(p, m) av * k / (p * dv + du) *
      F21(av * (1 + m) + 1, p * r + 1, p * r + 2, -k)
    (1 - th) * T(1, 0) + 2 * th * T(2, 0) + 2 * th * T(1, 1) - 4 * th * T(2, 1)
  }
  bd_cost <- function(e, av, dv, du, th) {              # Theorem 2, eps = u + v
    k <- k1f(e, av, dv); r <- dv / du
    T <- function(p, m) {
      A <- av * (1 + m)
      av * k^(-A) / (du * (A + p * r)) * F21(A + 1, A + p * r, A + p * r + 1, -1 / k)
    }
    (1 - th) * T(1, 0) + 2 * th * T(2, 0) + 2 * th * T(1, 1) - 4 * th * T(2, 1)
  }

  for (av in c(0.4, 1, 3)) {
    for (th in c(0, 0.7, -0.9)) {
      ep <- seq(-3, 2, length.out = 9)
      expect_equal(.cop_test_dens(ep, av, 0.4, 0.9, th, S = 1),
        bd_prod(ep, av, 0.4, 0.9, th), tolerance = 1e-10,
        info = paste("production", av, th))
      ec <- seq(-1, 4, length.out = 9)
      expect_equal(.cop_test_dens(ec, av, 0.4, 0.9, th, S = -1),
        bd_cost(ec, av, 0.4, 0.9, th), tolerance = 1e-10,
        info = paste("cost", av, th))
    }
  }
})

test_that("a cost frontier is not fitted with the density mirrored", {
  ## REGRESSION. Up to and including 1.2.0 the composition used eps + S*u where
  ## it needed eps + u. `eps` is already sign-normalized, so multiplying the
  ## quadrature node by S as well reflected the density: with inefdec = FALSE
  ## copsfm() maximized f_{v-u}(y - Xb) instead of f_{v+u}(y - Xb), and on
  ## clean cost data with sigma_u = 1, sigma_v = 0.4 it returned sigma_u =
  ## 0.031 and sigma_v = 0.69. The slopes survived; nothing else did.
  gl <- sfa:::.gauss_legendre_01(128L)
  tt <- gl$nodes; wt <- gl$weights; uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  su <- 1; sv <- 0.4
  x <- seq(-2, 3, length.out = 7)          # observed y - Xb, a COST frontier
  for (S in c(1, -1)) {
    eps <- S * x
    U <- matrix(uu, length(eps), length(uu), byrow = TRUE)
    vv <- S * (eps + U)
    lg <- dnorm(vv, 0, sv, log = TRUE) + log(2) + dnorm(U, 0, su, log = TRUE)
    got <- sfa:::.log_row_sum_exp(sweep(lg, 2, log(wt) + log(jac), "+"))
    ## The closed form, skewed in whichever direction the orientation implies:
    ## S = +1 is v - u, negatively skewed, so the tilt is Phi(-x lambda/sigma);
    ## S = -1 is v + u and the tilt reverses.
    sig <- sqrt(su^2 + sv^2); lam <- su / sv
    want <- log(2) - log(sig) + dnorm(x / sig, log = TRUE) +
      pnorm(-S * x * lam / sig, log.p = TRUE)
    expect_equal(got, want, tolerance = 1e-8, info = paste("S =", S))
  }
})

test_that("copsfm reports the parameters each (copula, udist, vdist) implies", {
  skip_on_cran()
  d <- cop_gen(4, 400, 0.3)
  nm <- function(...) rownames(copsfm(y ~ x1 + x2, data = d, n_nodes = 32, ...)$out)
  base <- c("(Intercept)", "x1", "x2", "sigma_u")
  ## "independent" estimates no dependence parameter at all, rather than
  ## carrying an unidentified one pinned at its independence value.
  expect_identical(nm(copula = "independent"), c(base, "sigma_v"))
  expect_identical(nm(copula = "fgm"), c(base, "sigma_v", "theta"))
  ## The generalized-logistic scale is delta_v, NOT sigma_v: it is not the
  ## standard deviation, and naming it sigma_v would invite reading it as one.
  expect_identical(nm(copula = "independent", vdist = "logistic"),
    c(base, "delta_v"))
  expect_identical(nm(copula = "independent", vdist = "glogistic"),
    c(base, "delta_v", "alpha_v"))
  expect_identical(nm(copula = "fgm", udist = "exponential", vdist = "glogistic"),
    c(base, "delta_v", "alpha_v", "theta"))

  f <- copsfm(y ~ x1 + x2, data = d, copula = "independent",
    udist = "exponential", vdist = "logistic", n_nodes = 32)
  expect_identical(f$udist, "exponential")
  expect_identical(f$vdist, "logistic")
  expect_true(is.na(f$copula_par))
  expect_equal(f$alpha_v, 1)                 # fixed, not estimated
  ## Battese-Coelli sits alongside the JLMS predictor, and is a probability.
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_equal(length(f$exp_u_hat), nrow(d))
  ## "independent" is a legitimate choice and must not warn about recovery.
  expect_silent(sfa:::.cop_warn_family("independent"))
})

test_that("the frontier slopes survive an exponential u and a skewed v", {
  skip_on_cran()
  ## FGM draws by the conditional method: with A = theta(1 - 2a), the
  ## conditional CDF b + A b(1 - b) is a quadratic with one root in [0, 1].
  set.seed(21)
  n <- 1500; th <- 0.6
  a <- runif(n); p <- runif(n); A <- th * (1 - 2 * a)
  b <- ifelse(abs(A) < 1e-10, p,
    ((1 + A) - sqrt(pmax((1 + A)^2 - 4 * A * p, 0))) / (2 * A))
  x1 <- rnorm(n); x2 <- rnorm(n)
  v <- sfa:::.gl_q(a, 2.5, 0.4)
  u <- qexp(b, rate = 1 / 0.9)
  d <- data.frame(y = 0.5 + 0.8 * x1 - 0.4 * x2 + v - u, x1 = x1, x2 = x2)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "fgm", udist = "exponential",
    vdist = "glogistic", n_nodes = 64)
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.08)
  expect_equal(unname(f$out["x2", "par"]), -0.4, tolerance = 0.08)
  expect_equal(unname(f$out["sigma_u", "par"]), 0.9, tolerance = 0.25)
})

test_that("the marginals compose correctly in the combinations the paper omits", {
  ## The Theorem 1/2 check above covers exponential u with generalized-logistic
  ## v, which is the pair Bonanno and Domma derive. The quadrature admits four
  ## pairs, and the other two are checked here against something independent --
  ## which also isolates the MARGINAL dispatch from the copula machinery, since
  ## both of these are at independence.
  gl <- sfa:::.gauss_legendre_01(256L)
  tt <- gl$nodes; wt <- gl$weights; uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  e <- seq(-4, 2, length.out = 11)
  su <- 0.9; sv <- 0.4
  U <- matrix(uu, length(e), length(uu), byrow = TRUE)
  compose <- function(lv, lu) {
    sfa:::.log_row_sum_exp(sweep(lv + lu, 2, log(wt) + log(jac), "+"))
  }

  ## normal v + exponential u IS the normal-exponential model, which has a
  ## closed form that sfm(model_name = "NE") uses.
  got <- compose(sfa:::.cop_v_ld(e + U, sv, 1, "normal"),
    sfa:::.cop_u_ld(U, su, "exponential"))
  want <- -log(su) + sv^2 / (2 * su^2) + e / su +
    pnorm(-e / sv - sv / su, log.p = TRUE)
  expect_equal(got, want, tolerance = 1e-12)

  ## half-normal u with a skewed noise has no closed form, so: integration.
  got2 <- compose(sfa:::.cop_v_ld(e + U, 0.4, 2.5, "glogistic"),
    sfa:::.cop_u_ld(U, su, "hnormal"))
  want2 <- vapply(e, function(ei) log(integrate(function(u)
    exp(sfa:::.gl_ld(ei + u, 2.5, 0.4)) * 2 * dnorm(u, 0, su), 0, Inf,
    rel.tol = 1e-12)$value), 0)
  expect_equal(got2, want2, tolerance = 1e-12)

  ## The CDFs the copula is fed are the marginals' own, not approximations.
  expect_equal(sfa:::.cop_u_p(c(0.1, 1, 5), 0.9, "exponential"),
    pexp(c(0.1, 1, 5), rate = 1 / 0.9))
  expect_equal(sfa:::.cop_u_p(c(0.1, 1, 5), 0.9, "hnormal"),
    2 * pnorm(c(0.1, 1, 5) / 0.9) - 1)
  expect_equal(sfa:::.cop_v_p(c(-2, 0, 3), 0.4, 1, "normal"),
    pnorm(c(-2, 0, 3) / 0.4))
})

test_that("a parameter on its bound does not cost the standard errors", {
  ## REGRESSION. copsfm()'s likelihood refused out-of-range draws with
  ## .Machine$double.xmax. optim() differences the objective to build its
  ## gradient, and differencing 1.8e308 overflows to a non-finite value, so the
  ## final stage aborted and every standard error came back NA. It bites
  ## whenever a parameter ENDS UP on a bound, which a dependence parameter
  ## routinely does -- FGM's theta reaches -0.99 on these data. Now a large
  ## FINITE penalty, and the two bounded parameters are clamped rather than
  ## refused, since the optimizers already hold them in range.
  skip_on_cran()
  d <- cop_gen(4, 400, 0.3)
  for (cfg in list(c("fgm", "exponential", "glogistic"),
                   c("gaussian", "exponential", "glogistic"),
                   c("frank", "exponential", "logistic"))) {
    f <- suppressWarnings(copsfm(y ~ x1 + x2, data = d, copula = cfg[1],
      udist = cfg[2], vdist = cfg[3], n_nodes = 32))
    expect_true(all(is.finite(f$out[, "st_err"])),
      info = paste(cfg, collapse = "/"))
  }
})
