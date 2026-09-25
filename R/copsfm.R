## Copula stochastic frontier: dependence between the noise and the
## inefficiency. Gap K5. See notes/code_history/copsfm.md.
##
## Every other model in this package assumes v and u are independent. That
## assumption is ubiquitous in the frontier literature and, as Smith (2008)
## argued, mostly an artefact of how easily the likelihood is then written
## down: if a farmer misjudges a seasonal shock and mis-allocates inputs in
## response, the shock and the inefficiency are the same event seen twice.
##
## The package's own taxonomy already states the general form (JSS paper, Eq 6):
##
##   f_{V,U}(v, u) = f_V(v) f_U(u) c(F_V(v), F_U(u); rho)
##
## with every existing specification setting c = 1. This file relaxes that, and
## nothing else: the marginals are the ordinary normal and half-normal.
##
## Composing, with eps = v - S*u so that v = eps + S*u,
##
##   f_eps(eps) = int_0^inf f_V(eps + S u) f_U(u) c(F_V(eps + S u), F_U(u)) du
##
## which is a one-dimensional integral per observation, done by Gauss-Legendre
## on a transformed range rather than by simulation.

## Copula densities on the unit square, in LOGS.
##
## Only families whose density can be written down exactly and checked are
## included. Each is verified in the tests two ways: it integrates to 1 over
## [0,1]^2, and at the independence parameter it returns exactly 1.
.cop_logc <- function(w1, w2, par, family) {
  switch(family,
    independent = rep(0, length(w1)),
    gaussian = {
      ## rho in (-1, 1). a = Phi^-1(w1), b = Phi^-1(w2):
      ##   c = (1-rho^2)^-1/2 exp{ [2 rho a b - rho^2 (a^2 + b^2)] / [2(1-rho^2)] }
      rho <- par
      if (abs(rho) >= 1) return(rep(NA_real_, length(w1)))
      a <- stats::qnorm(w1)
      b <- stats::qnorm(w2)
      om <- 1 - rho^2
      -0.5 * log(om) + (2 * rho * a * b - rho^2 * (a^2 + b^2)) / (2 * om)
    },
    fgm = {
      ## Farlie-Gumbel-Morgenstern: c = 1 + theta (1-2 w1)(1-2 w2), |theta| <= 1.
      ## Bounded dependence (Spearman rho is theta/3), so it cannot represent
      ## strong association -- which is exactly why it is a useful contrast to
      ## the Gaussian rather than a redundant second option.
      th <- par
      if (abs(th) > 1) return(rep(NA_real_, length(w1)))
      v <- 1 + th * (1 - 2 * w1) * (1 - 2 * w2)
      ifelse(v > 0, log(v), NA_real_)
    },
    ## Frank: the only other one-parameter family that spans BOTH signs of
    ## dependence over the full range. theta and (1 - e^-theta) always share a
    ## sign, so their product is positive and can be logged directly; the
    ## independence limit theta -> 0 is taken explicitly because the ratio is
    ## 0/0 there.
    frank = {
      th <- par
      if (!is.finite(th)) return(rep(NA_real_, length(w1)))
      if (abs(th) < 1e-8) return(rep(0, length(w1)))
      em <- -expm1(-th)                      # 1 - exp(-th), signed like th
      D <- em - (-expm1(-th * w1)) * (-expm1(-th * w2))
      ok <- is.finite(D) & D != 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- log(th * em) - th * (w1[ok] + w2[ok]) - 2 * log(abs(D[ok]))
      out
    },
    ## Clayton: lower-tail dependence, theta > 0, independence as theta -> 0.
    clayton = {
      th <- par
      if (!is.finite(th) || th <= 0) {
        return(if (isTRUE(all.equal(th, 0))) rep(0, length(w1)) else rep(NA_real_, length(w1)))
      }
      z <- w1^(-th) + w2^(-th) - 1
      ok <- is.finite(z) & z > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- log1p(th) - (1 + th) * (log(w1[ok]) + log(w2[ok])) -
        (2 + 1 / th) * log(z[ok])
      out
    },
    ## Gumbel: upper-tail dependence, theta >= 1, independence at theta = 1.
    gumbel = {
      th <- par
      if (!is.finite(th) || th < 1) return(rep(NA_real_, length(w1)))
      if (abs(th - 1) < 1e-12) return(rep(0, length(w1)))
      x <- -log(w1); y <- -log(w2)
      A <- (x^th + y^th)^(1 / th)
      ok <- is.finite(A) & A > 0 & x > 0 & y > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- -A[ok] + (th - 1) * (log(x[ok]) + log(y[ok])) +
        (1 - 2 * th) * log(A[ok]) + log(A[ok] + th - 1) -
        log(w1[ok]) - log(w2[ok])
      out
    },
    ## Joe: upper-tail dependence, heavier than Gumbel; theta >= 1.
    joe = {
      th <- par
      if (!is.finite(th) || th < 1) return(rep(NA_real_, length(w1)))
      if (abs(th - 1) < 1e-12) return(rep(0, length(w1)))
      a <- (1 - w1)^th; b <- (1 - w2)^th
      z <- a + b - a * b
      ok <- is.finite(z) & z > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- (1 / th - 2) * log(z[ok]) +
        (th - 1) * (log1p(-w1[ok]) + log1p(-w2[ok])) +
        log(th - 1 + z[ok])
      out
    },
    NA_real_
  )
}

## Rotations. Clayton, Gumbel and Joe carry only POSITIVE dependence, which is
## a real restriction here: nothing rules out a negative association between
## noise and inefficiency, and a family that cannot express one will report the
## independence boundary instead of a negative estimate. The 90 and 270 degree
## rotations are exact reflections of the density, so they need no new algebra:
##   c_90(u,v)  = c(1-u, v)      c_180(u,v) = c(1-u, 1-v)      c_270 = c(u, 1-v)
.COP_ROT <- c(gaussian = NA, fgm = NA, frank = NA, independent = NA,
  clayton = 0, gumbel = 0, joe = 0,
  clayton90 = 90, clayton180 = 180, clayton270 = 270,
  gumbel90 = 90, gumbel180 = 180, gumbel270 = 270,
  joe90 = 90, joe180 = 180, joe270 = 270)

.cop_base <- function(family) sub("(90|180|270)$", "", family)

.cop_logc_rot <- function(w1, w2, par, family) {
  rot <- .COP_ROT[[family]]
  base <- .cop_base(family)
  if (is.null(rot) || is.na(rot) || rot == 0) return(.cop_logc(w1, w2, par, base))
  if (rot == 90) return(.cop_logc(1 - w1, w2, par, base))
  if (rot == 180) return(.cop_logc(1 - w1, 1 - w2, par, base))
  .cop_logc(w1, 1 - w2, par, base)
}

## Independence value of each family's parameter, and its admissible range.
.cop_spec <- function(family) {
  ## Upper limits are where the density stops being computable in double
  ## precision, not where the family stops being defined: Clayton at theta = 28
  ## already has a Kendall tau of 0.93, and past that u^-theta overflows.
  ## Independence sits ON the lower bound for Clayton/Gumbel/Joe, so par0 is set
  ## just inside it -- starting the optimiser exactly at a bound is how a search
  ## reports the boundary back as an estimate.
  base <- .cop_base(family)
  switch(base,
    gaussian = list(par0 = 0, lo = -0.95, hi = 0.95, name = "rho"),
    fgm      = list(par0 = 0, lo = -0.99, hi = 0.99, name = "theta"),
    frank    = list(par0 = 1e-4, lo = -35, hi = 35, name = "theta"),
    clayton  = list(par0 = 0.05, lo = 1e-6, hi = 28, name = "theta"),
    gumbel   = list(par0 = 1.0001, lo = 1, hi = 17, name = "theta"),
    joe      = list(par0 = 1.0001, lo = 1, hi = 30, name = "theta"),
    list(par0 = 0, lo = -0.95, hi = 0.95, name = "par")
  )
}

## Which families have been shown to recover their own dependence parameter.
##
## Measured 2026-09-04: 25 samples per family generated FROM that family at
## n = 400, refitted with it, counting how often the estimate came back on the
## independence boundary. A correct density does not imply a recoverable
## parameter, and for most of these families it is not recoverable: on Gumbel
## data with Spearman 0.685 at n = 2000, every family INCLUDING the true one
## returns the boundary and their log-likelihoods differ by less than 0.04.
## The likelihood is flat in the dependence parameter. That is the model, not
## the code -- each density is verified against the second mixed partial of its
## own CDF in test-copsfm.R.
.COP_COLLAPSE <- c(gumbel = 36, joe = 40, clayton270 = 56, gumbel90 = 60)
.COP_VERIFIED <- c("gaussian", "fgm", "frank", "clayton")

.cop_warn_family <- function(family) {
  ## "independent" has no parameter to recover, so there is nothing to warn about.
  if (identical(family, "independent") || family %in% .COP_VERIFIED) {
    return(invisible(NULL))
  }
  ## `[` not `[[`: a name absent from the vector must give NA, not an error.
  rate <- unname(.COP_COLLAPSE[family])
  if (!is.na(rate)) {
    warning("copsfm(copula = \"", family, "\"): this family did not recover its ",
      "own dependence parameter in testing -- on 25 samples generated from it ",
      "at n = 400, ", rate, "% of fits returned the independence boundary. The ",
      "density is correct; the parameter is weakly identified, and the ",
      "likelihood is close to flat in it. Prefer copula = \"frank\" or ",
      "\"clayton\", which recovered with no boundary collapses, and treat any ",
      "estimate from this family as exploratory. See ?copsfm.",
      call. = FALSE)
  } else {
    warning("copsfm(copula = \"", family, "\"): this family's density is ",
      "verified but its RECOVERY has not been measured. The families that were ",
      "measured collapsed to the independence boundary on 36-60% of samples ",
      "generated from themselves, so do not assume this one behaves better. ",
      "Prefer copula = \"frank\" or \"clayton\". See ?copsfm.",
      call. = FALSE)
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Marginal distributions for the two error components. Gap L18.
## ---------------------------------------------------------------------------
##
## Bonanno and Domma (2022) pair an EXPONENTIAL inefficiency with a GENERALIZED
## LOGISTIC noise carrying its own skewness parameter. Their argument is that
## the wrong-skewness anomaly is a specification failure rather than a small
## sample accident: the third moment of the composed error (their Eq. 5)
## depends on the skew of v and on the dependence between v and u, and not only
## on the skew of u, so a model in which v cannot be skewed has nowhere to put
## a positive residual skewness except in a sigma_u that collapses to zero.
##
## Nothing here is copula-specific. The convolution is done by quadrature, so
## the marginals vary independently of the dependence structure, and the
## paper's four specifications -- (I,S), (I,A), (D,S), (D,A) -- are four
## (copula, vdist) pairs rather than four separate estimators.

## Generalized logistic GL(alpha, delta), in Bonanno and Domma's centred
## parameterization:
##
##   z    = v / delta + [psi(alpha) - psi(1)]
##   G(v) = (1 + e^-z)^-alpha
##   g(v) = (alpha / delta) e^-z (1 + e^-z)^-(alpha+1)
##
## The location term is the whole point: it holds E[V] = 0 for EVERY (alpha,
## delta), so the two parameters move the spread and the skew of the noise
## without moving its mean. A noise term with a non-zero mean is not noise, it
## is part of the frontier, and would not be separately identified from the
## intercept.
##
## Then Var(V) = delta^2 [psi'(alpha) + psi'(1)] -- NOT delta^2, which is why
## the fitted scale is reported as `delta_v` and not as `sigma_v` -- and
## E[V - E V]^3 = delta^3 [psi''(alpha) - psi''(1)]. alpha = 1 is the ordinary
## logistic with scale delta; alpha < 1 skews negative, alpha > 1 positive.
.gl_z <- function(v, alpha, delta) v / delta + (digamma(alpha) - digamma(1))

.gl_ld <- function(v, alpha, delta) {
  z <- .gl_z(v, alpha, delta)
  log(alpha) - log(delta) - z - (alpha + 1) * .log1pexp(-z)
}

.gl_lp <- function(v, alpha, delta) -alpha * .log1pexp(-.gl_z(v, alpha, delta))

.gl_q <- function(p, alpha, delta) {
  delta * (-log(p^(-1 / alpha) - 1) - (digamma(alpha) - digamma(1)))
}

## alpha is bounded away from 0 rather than merely from below. psi(alpha) ~
## -1/alpha as alpha -> 0, so the centring shift is delta/alpha and the density
## walks off to infinity while remaining perfectly finite at every point the
## optimizer evaluates -- a silent failure, not an error.
.GL_ALPHA <- c(0.05, 20)

## (log density, CDF) of each marginal. Both take the SCALE parameter; the
## shape argument is ignored except by the generalized logistic.
.cop_u_ld <- function(u, su, udist) {
  if (identical(udist, "exponential")) {
    stats::dexp(u, rate = 1 / su, log = TRUE)
  } else {
    log(2) + stats::dnorm(u, 0, su, log = TRUE)
  }
}

.cop_u_p <- function(u, su, udist) {
  if (identical(udist, "exponential")) -expm1(-u / su) else 2 * stats::pnorm(u / su) - 1
}

.cop_v_ld <- function(v, sv, av, vdist) {
  if (identical(vdist, "normal")) stats::dnorm(v, 0, sv, log = TRUE) else .gl_ld(v, av, sv)
}

.cop_v_p <- function(v, sv, av, vdist) {
  if (identical(vdist, "normal")) stats::pnorm(v / sv) else exp(.gl_lp(v, av, sv))
}

## Moments the starting values need: the coefficient of sigma_u^3 in the third
## central moment of -u, the coefficient of sigma_u^2 in Var(u), E[u]/sigma_u,
## and the coefficient of delta_v^2 in Var(v) at alpha = 1.
.cop_u_k3 <- function(udist) if (identical(udist, "exponential")) -2 else sqrt(2 / pi) * (1 - 4 / pi)
.cop_u_v2 <- function(udist) if (identical(udist, "exponential")) 1 else 1 - 2 / pi
.cop_u_m1 <- function(udist) if (identical(udist, "exponential")) 1 else sqrt(2 / pi)
.cop_v_v2 <- function(vdist) if (identical(vdist, "normal")) 1 else pi^2 / 3

copsfm <- function(formula,
                   data,
                   copula = c("gaussian", "fgm", "frank",
                              "clayton", "clayton90", "clayton180", "clayton270",
                              "gumbel", "gumbel90", "gumbel180", "gumbel270",
                              "joe", "joe90", "joe180", "joe270",
                              "independent"),
                   inefdec = TRUE,
                   udist = c("hnormal", "exponential"),
                   vdist = c("normal", "logistic", "glogistic"),
                   n_nodes = 128,
                   maxit.bobyqa = 10000,
                   maxit.psoptim = 1000,
                   maxit.optim = 1000,
                   start_val = FALSE,
                   PSopt = FALSE,
                   optHessian = TRUE,
                   Method = "L-BFGS-B",
                   verbose = FALSE,
                   keep_objective = FALSE,
                   rand.psoptim = NULL) {
  call <- match.call()
  copula <- match.arg(copula)
  udist <- match.arg(udist)
  vdist <- match.arg(vdist)
  .cop_warn_family(copula)
  has_cop <- !identical(copula, "independent")
  Start.Time <- Sys.time()
  cz <- .SFA_CONSTANTS

  if (missing(formula) || !inherits(formula, "formula") || length(formula) != 3L) {
    stop("copsfm(): `formula` must be a two-sided formula, e.g. y ~ x1 + x2.",
      call. = FALSE
    )
  }
  if (any(grepl("|", deparse(formula), fixed = TRUE))) {
    stop("copsfm(): `formula` must not contain a `|` segment. Determinants of ",
      "the variance are not yet supported alongside a copula.",
      call. = FALSE
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("copsfm(): `data` must be a data.frame.", call. = FALSE)
  }
  if (length(inefdec) != 1L || !is.logical(inefdec) || is.na(inefdec)) {
    stop("copsfm(): `inefdec` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(n_nodes) != 1L || !is.numeric(n_nodes) || n_nodes < 16) {
    stop("copsfm(): `n_nodes` must be a single number >= 16.", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(stats::terms(mf), mf)
  n <- length(y)
  k <- ncol(X)
  np <- k + 2L + (vdist == "glogistic") + has_cop
  if (n <= np) {
    stop("copsfm(): ", n, " observations cannot identify ", np,
      " parameters.",
      call. = FALSE
    )
  }
  S <- if (isTRUE(inefdec)) 1 else -1
  cs <- if (has_cop) .cop_spec(copula) else NULL

  ## Gauss-Legendre nodes on (0, 1), mapped to u in (0, Inf) by u = t/(1-t).
  ## The map is chosen so the integrand's mass -- which sits at moderate u --
  ## lands in the middle of the node range rather than in a tail.
  ##
  ## 128 nodes, not 64. Measured against the closed-form normal/half-normal
  ## density at independence, the maximum absolute error in log f is 1.3e-2 at
  ## 32 nodes, 2.5e-5 at 64 and 8.2e-14 at 128. 64 looks adequate and is not:
  ## an error of 1e-5 per observation is an error of 0.02 in a log-likelihood
  ## over 2000 points, which is the scale at which modes are being compared.
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  tt <- gl$nodes
  wt <- gl$weights
  uu <- tt / (1 - tt)
  jac <- 1 / (1 - tt)^2

  ## The per-node log joint density, rows observations and columns quadrature
  ## nodes. Shared by the likelihood and by the efficiency predictors, so the
  ## objective and the prediction cannot come to disagree about the model.
  ##
  ## `eps` is S * (y - Xb), which is distributed as (S v) - u. Adding back the
  ## node therefore recovers S*v, and the NOISE ITSELF is S times that. That
  ## last multiplication is not cosmetic: for a cost frontier the sign-
  ## normalized error is (-v) - u, and a SKEWED v does not survive the
  ## reflection the way a normal one does. The same reflection is why F_V is
  ## evaluated at the actual v -- so that a positive dependence parameter means
  ## the same thing, association between the real v and the real u, in both
  ## orientations.
  .node_lg <- function(eps, su, sv, av, cpar, U) {
    vv <- S * (eps + U)
    lg <- .cop_v_ld(vv, sv, av, vdist) + .cop_u_ld(U, su, udist)
    if (has_cop && !is.null(cpar)) {
      ## Marginal CDFs, clamped off the endpoints: qnorm(0) is -Inf.
      w1 <- pmin(pmax(.cop_v_p(vv, sv, av, vdist), 1e-12), 1 - 1e-12)
      w2 <- pmin(pmax(.cop_u_p(U, su, udist), 1e-12), 1 - 1e-12)
      lc <- .cop_logc_rot(as.numeric(w1), as.numeric(w2), cpar, copula)
      if (any(!is.finite(lc))) return(NULL)
      lg <- lg + matrix(lc, nrow = nrow(lg))
    }
    sweep(lg, 2, log(wt) + log(jac), "+")
  }

  ## log f_eps(eps_i), the composed density, by quadrature over u.
  .log_dens <- function(eps, su, sv, av, cpar) {
    U <- matrix(uu, nrow = length(eps), ncol = length(uu), byrow = TRUE)
    lg <- .node_lg(eps, su, sv, av, cpar, U)
    if (is.null(lg)) return(rep(NA_real_, length(eps)))
    .log_row_sum_exp(lg)
  }

  ## Starting values: OLS plus Olson's moments, independence for the copula and
  ## a symmetric noise. The moment identities are per-distribution -- Olson
  ## reads sigma_u out of the residual skewness, and the exponential and the
  ## half-normal put different multiples of sigma_u^3 there.
  ols <- stats::lm.fit(x = X, y = y)
  e0 <- as.numeric(ols$residuals)
  m2 <- mean(e0^2); m3 <- mean(e0^3)
  k3 <- .cop_u_k3(udist)
  su0 <- if (S * m3 < 0) (S * m3 / k3)^(1 / 3) else 0.5 * stats::sd(e0)
  su0 <- if (is.finite(su0)) max(su0, 1e-3) else max(0.5 * stats::sd(e0), 1e-3)
  ## Var(eps) = Var(v) + Var(u), and for a generalized logistic at alpha = 1
  ## Var(v) is delta_v^2 pi^2/3 rather than delta_v^2.
  sv0 <- sqrt(max(m2 - .cop_u_v2(udist) * su0^2, 1e-6) / .cop_v_v2(vdist))
  b0 <- as.numeric(ols$coefficients)
  if (attr(stats::terms(mf), "intercept") == 1L) {
    b0[1L] <- b0[1L] + S * su0 * .cop_u_m1(udist)
  }

  ## Layout: beta, log sigma_u, log delta_v, then log alpha_v where the noise
  ## carries a shape, then the copula parameter where there is a copula.
  i_b <- seq_len(k); i_su <- k + 1L; i_sv <- k + 2L
  i_av <- NA_integer_; i_c <- NA_integer_
  start_v <- c(b0, log(su0), log(sv0))
  par_names <- c(colnames(X), "sigma_u",
    if (identical(vdist, "normal")) "sigma_v" else "delta_v")
  span <- pmax(10 * abs(b0), 10)
  lower1 <- c(b0 - span, log(1e-6), log(1e-6))
  upper1 <- c(b0 + span, log(1e4), log(1e4))
  if (identical(vdist, "glogistic")) {
    i_av <- length(start_v) + 1L
    start_v <- c(start_v, 0)                    # alpha_v = 1, the symmetric case
    par_names <- c(par_names, "alpha_v")
    lower1 <- c(lower1, log(.GL_ALPHA[1L])); upper1 <- c(upper1, log(.GL_ALPHA[2L]))
  }
  if (has_cop) {
    i_c <- length(start_v) + 1L
    start_v <- c(start_v, cs$par0)
    par_names <- c(par_names, cs$name)
    lower1 <- c(lower1, cs$lo); upper1 <- c(upper1, cs$hi)
  }
  if (isTRUE(start_val)) names(start_v) <- par_names

  ## A large FINITE penalty, not .Machine$double.xmax. optim() differences the
  ## objective to build its gradient, and differencing 1.8e308 overflows to a
  ## non-finite value -- which aborts the final stage with "non-finite
  ## finite-difference value" and costs the standard errors for the whole fit.
  ## sfm()'s NGE, NLN and NW branches already use 1e12 for exactly this reason.
  .PEN <- 1e12

  ## `per_obs = TRUE` returns the vector of per-observation log-likelihood
  ## contributions instead of the negative sum, for estfun.sfareg() and
  ## vcov(type = "bhhh"). The optimizers call this with one argument.
  ## Note the scale: th holds LOG sigma_u and LOG sigma_v, so the scores are
  ## on the log scale and vcov() delta-corrects them via `par_scale` below.
  like.fn <- function(th, per_obs = FALSE) {
    ## A refused draw must keep its length under per_obs; the scalar barrier
    ## would collapse estfun()'s matrix to a single row.
    .bail <- function() if (isTRUE(per_obs)) rep(-.PEN / n, n) else .PEN
    if (!all(is.finite(th))) return(.bail())
    su <- exp(pmin(th[i_su], 12)); sv <- exp(pmin(th[i_sv], 12))
    ## CLAMPED, not rejected. The optimizers already hold th[i_av] inside these
    ## bounds; returning a huge penalty for the last ulp of floating-point slop
    ## at the bound is how the final optim() stage aborts with a non-finite
    ## finite-difference value, which is exactly what it did.
    av <- if (is.na(i_av)) 1 else {
      exp(min(max(th[i_av], log(.GL_ALPHA[1L])), log(.GL_ALPHA[2L])))
    }
    ## CLAMPED, not rejected -- the same reasoning as alpha_v above, and it
    ## matters more here because a dependence parameter routinely ENDS UP on
    ## its bound (FGM's theta especially). A 1e12 cliff one finite-difference
    ## step from the reported optimum costs the standard errors.
    cpar <- if (has_cop) min(max(th[i_c], cs$lo), cs$hi) else NULL
    eps <- S * as.numeric(y - X %*% th[i_b])
    ll <- .log_dens(eps, su, sv, av, cpar)
    if (any(!is.finite(ll))) return(.bail())
    if (isTRUE(per_obs)) return(ll)
    ## The scaffold MINIMIZES; every likelihood here returns the negative sum.
    -sum(ll)
  }

  Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v,
    lower.bobyqa = lower1, upper.bobyqa = upper1,
    maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, rhobeg = NA, rhoend = NA,
    verbose = verbose
  )
  start_v <- Opt.Bobyqa$start_v; bob1 <- Opt.Bobyqa$bob1

  Opt.Psoptim <- opt.psoptim(fn = like.fn, start_v,
    lower.psoptim = lower1, rand.psoptim = rand.psoptim,
    upper.psoptim = upper1, maxit.psoptim, psopt.TF = PSopt,
    rand.order = FALSE, verbose = verbose
  )
  start_v <- Opt.Psoptim$start_v; opt00 <- Opt.Psoptim$opt00

  Opt.Optim <- opt.optim(fn = like.fn, start_v = start_v,
    lower.optim = lower1, upper.optim = upper1, maxit.optim = maxit.optim,
    opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
  )
  start_v <- Opt.Optim$start_v; opt <- Opt.Optim$opt
  End.Time <- end.time(Start.Time)

  if (optHessian == FALSE & PSopt == FALSE) { opt <- bob1; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == FALSE & PSopt == TRUE)  { opt <- opt00; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == TRUE) {
    st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(diag(solve(opt$hessian))))
    }
  }

  th <- opt$par
  su <- exp(th[i_su]); sv <- exp(th[i_sv])
  av <- if (is.na(i_av)) 1 else exp(th[i_av])
  cpar <- if (has_cop) th[i_c] else NA_real_
  par <- c(th[i_b], su, sv)
  se <- c(st_err[i_b], su * st_err[i_su], sv * st_err[i_sv])
  ## alpha_v is estimated as log alpha_v and reported on its own scale, so its
  ## standard error carries the same delta-method factor the two scales do.
  if (!is.na(i_av)) { par <- c(par, av); se <- c(se, av * st_err[i_av]) }
  if (has_cop) { par <- c(par, cpar); se <- c(se, st_err[i_c]) }

  out <- matrix(NA_real_, 3L, length(par))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- par_names
  out[1, ] <- par; out[2, ] <- se; out[3, ] <- par / se

  ## The efficiency predictors, by the same quadrature -- and from the same
  ## closure -- that the likelihood used.
  eps <- S * as.numeric(y - X %*% th[i_b])
  Um <- matrix(uu, nrow = n, ncol = length(uu), byrow = TRUE)
  lg <- .node_lg(eps, su, sv, av, cpar, Um)
  wts <- exp(lg - .log_row_sum_exp(lg))
  jlms <- as.numeric(rowSums(wts * Um))
  ## Battese-Coelli, E[exp(-u) | eps]. Free here, and it is the quantity the
  ## copula literature reports; exp(-jlms) is kept as `efficiency` because that
  ## is what this function has always returned under that name.
  exp_u_hat <- as.numeric(rowSums(wts * exp(-Um)))

  results <- list(
    t(out), c(opt), End.Time, start_v, "COP", formula, copula, cpar,
    udist, vdist, av, jlms, exp(-jlms), exp_u_hat, S, n, as.integer(n_nodes),
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c(
    "out", "opt", "total_time", "start_v", "model_name", "formula", "copula",
    "copula_par", "udist", "vdist", "alpha_v", "jlms", "efficiency",
    "exp_u_hat", "S", "nobs", "n_nodes",
    "coefficients", "std.errors", "t.values", "call"
  )
  ## Optionally retain the objective for estfun()/vcov(type = "bhhh"), with
  ## the map from the estimation-scale vector `th` onto the reported one.
  ## These mirror `par` and `se` above exactly: same positions, and the same
  ## delta-method factors the Hessian standard errors already carry.
  results$par_index <- c(i_b, i_su, i_sv,
    if (!is.na(i_av)) i_av else NULL,
    if (has_cop) i_c else NULL
  )
  results$par_scale <- c(rep(1, length(i_b)), su, sv,
    if (!is.na(i_av)) av else NULL,
    if (has_cop) 1 else NULL
  )
  if (isTRUE(keep_objective)) results$objective <- like.fn
  results
}
