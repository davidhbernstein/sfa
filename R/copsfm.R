## Copula stochastic frontier: dependence between the noise and the
## inefficiency. Gap K5. See notes/code_history/copsfm.md.
##
## Every other model in this package assumes v and u are independent. That
## assumption is ubiquitous in the frontier literature and, as Smith (2008)
## argued, mostly an artefact of how easily the likelihood is then written
## down: if a farmer misjudges a seasonal shock and mis-allocates inputs in
## response, the shock and the inefficiency are the same event seen twice.
##
## The package's own taxonomy already states the general form (companion paper, Eq 6):
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
## q1, q2 are 1 - w1, 1 - w2 computed without cancellation; the upper-tail
## families need them, since 1 - w is only good to 1e-16 absolute (A74).
.cop_logc <- function(w1, w2, par, family, q1 = 1 - w1, q2 = 1 - w2) {
  switch(family,
    independent = rep(0, length(w1)),
    gaussian = {
      ## rho in (-1, 1). a = Phi^-1(w1), b = Phi^-1(w2):
      ##   c = (1-rho^2)^-1/2 exp{ [2 rho a b - rho^2 (a^2 + b^2)] / [2(1-rho^2)] }
      rho <- par
      if (abs(rho) >= 1) return(rep(NA_real_, length(w1)))
      a <- stats::qnorm(pmin(w1, q1)); a[w1 >= 0.5] <- -a[w1 >= 0.5]
      b <- stats::qnorm(pmin(w2, q2)); b[w2 >= 0.5] <- -b[w2 >= 0.5]
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
      l1 <- log(w1); l2 <- log(w2)
      z <- exp(-th * l1) + exp(-th * l2) - 1
      out <- log1p(th) - (1 + th) * (l1 + l2) - (2 + 1 / th) * log(z)
      out[!(is.finite(z) & z > 0)] <- NA_real_
      out
    },
    ## Gumbel: upper-tail dependence, theta >= 1, independence at theta = 1.
    gumbel = {
      th <- par
      if (!is.finite(th) || th < 1) return(rep(NA_real_, length(w1)))
      if (abs(th - 1) < 1e-12) return(rep(0, length(w1)))
      x <- -log(w1); h <- w1 >= 0.5; x[h] <- -log1p(-q1[h])
      y <- -log(w2); h <- w2 >= 0.5; y[h] <- -log1p(-q2[h])
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
      a <- q1^th; b <- q2^th
      z <- a + b - a * b
      ok <- is.finite(z) & z > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- (1 / th - 2) * log(z[ok]) +
        (th - 1) * (log(q1[ok]) + log(q2[ok])) +
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

.cop_logc_rot <- function(w1, w2, par, family, q1 = 1 - w1, q2 = 1 - w2) {
  rot <- .COP_ROT[[family]]
  base <- .cop_base(family)
  if (is.null(rot) || is.na(rot) || rot == 0) return(.cop_logc(w1, w2, par, base, q1, q2))
  if (rot == 90) return(.cop_logc(q1, w2, par, base, w1, q2))
  if (rot == 180) return(.cop_logc(q1, q2, par, base, w1, w2))
  .cop_logc(w1, q2, par, base, q1, w2)
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

## (F_U, 1 - F_U), each computed without cancellation. 2 Phi(x) - 1 keeps
## only absolute accuracy near 0, where a lower-tail copula is most sensitive,
## so below x = 0.25 the erf series is used; and 1 - F is never formed by
## subtraction, because the upper-tail families need it to relative accuracy
## and 1 - w is good only to 1e-16 absolute (A74).
.cop_u_pq <- function(u, su, udist) {
  x <- u / su
  if (identical(udist, "exponential")) return(list(p = -expm1(-x), q = exp(-x)))
  q <- 2 * stats::pnorm(-x)
  p <- 1 - q
  sm <- which(x < 0.25)
  if (length(sm)) {
    xs <- x[sm]; t <- xs * xs / 2
    p[sm] <- sqrt(2 / pi) * xs * (1 + t * (-1 / 3 + t * (1 / 10 + t * (-1 / 42 +
      t * (1 / 216 + t * (-1 / 1320 + t * (1 / 9360 + t * (-1 / 75600))))))))
  }
  list(p = p, q = q)
}
.cop_u_p <- function(u, su, udist) .cop_u_pq(u, su, udist)$p
.cop_u_q <- function(u, su, udist) .cop_u_pq(u, su, udist)$q

.cop_v_ld <- function(v, sv, av, vdist) {
  switch(vdist,
    normal = stats::dnorm(v, 0, sv, log = TRUE),
    laplace = -log(2 * sv) - abs(v) / sv,
    .gl_ld(v, av, sv)
  )
}

## (F_V, 1 - F_V), likewise: the smaller tail is computed, the other is its
## complement.
.cop_v_pq <- function(v, sv, av, vdist) {
  if (!(vdist %in% c("normal", "laplace"))) {
    lp <- .gl_lp(v, av, sv)
    return(list(p = exp(lp), q = -expm1(lp)))
  }
  r <- if (identical(vdist, "normal")) stats::pnorm(-abs(v) / sv) else 0.5 * exp(-abs(v) / sv)
  p <- r; q <- 1 - r
  up <- which(v > 0)
  p[up] <- q[up]; q[up] <- r[up]
  list(p = p, q = q)
}
.cop_v_p <- function(v, sv, av, vdist) .cop_v_pq(v, sv, av, vdist)$p
.cop_v_q <- function(v, sv, av, vdist) .cop_v_pq(v, sv, av, vdist)$q

## Moments the starting values need: the coefficient of sigma_u^3 in the third
## central moment of -u, the coefficient of sigma_u^2 in Var(u), E[u]/sigma_u,
## and the coefficient of delta_v^2 in Var(v) at alpha = 1.
.cop_u_k3 <- function(udist) if (identical(udist, "exponential")) -2 else sqrt(2 / pi) * (1 - 4 / pi)
.cop_u_v2 <- function(udist) if (identical(udist, "exponential")) 1 else 1 - 2 / pi
.cop_u_m1 <- function(udist) if (identical(udist, "exponential")) 1 else sqrt(2 / pi)
.cop_v_v2 <- function(vdist) switch(vdist, normal = 1, laplace = 2, pi^2 / 3)

## ---------------------------------------------------------------------------
## The integral over u (A74). See notes/code_history/copsfm.md.
## ---------------------------------------------------------------------------

## Log integrand of f_eps at nodes U, one eps per node. NULL when the copula
## density is not computable at the parameter.
.cop_node_ld <- function(U, eps, S, su, sv, av, udist, vdist, copula, cpar,
                         parts = FALSE) {
  vv <- S * (eps + U)
  lv <- .cop_v_ld(vv, sv, av, vdist)
  lu <- .cop_u_ld(U, su, udist)
  if (identical(copula, "independent") || is.null(cpar) || is.na(cpar)) {
    return(if (parts) list(lv = lv, lu = lu, lc = 0 * lu) else lv + lu)
  }
  cl <- function(w) { w[w < 1e-12] <- 1e-12; w[w > 1 - 1e-12] <- 1 - 1e-12; w }
  a <- .cop_v_pq(vv, sv, av, vdist)
  b <- .cop_u_pq(U, su, udist)
  lc <- .cop_logc_rot(cl(a$p), cl(b$p), cpar, copula, cl(a$q), cl(b$q))
  if (any(!is.finite(lc))) return(NULL)
  if (parts) list(lv = lv, lu = lu, lc = lc) else lv + lu + lc
}

## Log-odds of the marginals, used only to place breakpoints.
.cop_u_lo <- function(u, su, udist) {
  x <- u / su
  lq <- if (identical(udist, "exponential")) -x else log(2) + stats::pnorm(-x, log.p = TRUE)
  log(.cop_u_p(u, su, udist)) - lq
}
.cop_v_lo <- function(v, sv, av, vdist) {
  switch(vdist,
    normal = stats::pnorm(v / sv, log.p = TRUE) - stats::pnorm(-v / sv, log.p = TRUE),
    laplace = {
      x <- -abs(v) / sv
      sign(v) * (log1p(-0.5 * exp(x)) - log(0.5) - x)
    },
    {
      lp <- .gl_lp(v, av, sv)
      lp - log(-expm1(lp))
    }
  )
}

## Breakpoint levels, in log-odds. Dense where a marginal's bulk is, sparse in
## its tails, and +-27.63 is where the copula arguments are clamped at 1e-12.
.COP_LEV <- local({
  k <- log(1e-12) - log1p(-1e-12)
  l <- c(1, 2.5, 4.5, 7, 10, 14, 19, 24, -k, 34)
  c(rev(-l), 0, l)
})

.COP_REL <- c(1, 2.5, 4.5, 7, 10, 14, 19, 24, 34)

## u (resp. v) at which the log-odds of F_U (resp. F_V) equal each level.
.cop_u_at_lev <- function(su, udist) {
  a <- .COP_LEV
  if (identical(udist, "exponential")) return(su * ifelse(a > 30, a, log1p(exp(a))))
  lp <- stats::plogis(a, log.p = TRUE); lq <- stats::plogis(-a, log.p = TRUE)
  su * sqrt(ifelse(a < 0, stats::qchisq(lp, 1, log.p = TRUE),
    stats::qchisq(lq, 1, lower.tail = FALSE, log.p = TRUE)))
}
.cop_v_at_lev <- function(sv, av, vdist, b = .COP_LEV) {
  lp <- stats::plogis(b, log.p = TRUE); lq <- stats::plogis(-b, log.p = TRUE)
  switch(vdist,
    normal = sv * ifelse(b < 0, stats::qnorm(lp, log.p = TRUE), -stats::qnorm(lq, log.p = TRUE)),
    laplace = sv * ifelse(b < 0, log(2) + lp, ifelse(b == 0, 0, -(log(2) + lq))),
    {
      x <- -lp / av
      sv * (-ifelse(x > 30, x + log1p(-exp(-x)), log(expm1(x))) - (digamma(av) - digamma(1)))
    }
  )
}

## How fast the log copula density can change per unit of log-odds of either
## argument; it sets how finely a piece is cut.
.cop_kappa <- function(copula, cpar) {
  if (identical(copula, "independent") || is.null(cpar) || is.na(cpar)) return(0)
  switch(.cop_base(copula),
    gaussian = 2 * abs(cpar) / (1 - cpar^2),
    fgm = 1,
    frank = abs(cpar) / 2,
    clayton = 1 + cpar,
    cpar
  )
}

## Gauss-Legendre on (0, 1) for any n, including the 2- and 3-point rules
## .gauss_legendre_01() does not serve.
.cop_gl01 <- function(n) {
  if (n >= 5L) return(.gauss_legendre_01(as.integer(n)))
  x <- switch(n, 0, c(-1, 1) / sqrt(3), c(-1, 0, 1) * sqrt(3 / 5),
    c(-1, -1, 1, 1) * sqrt(3 / 7 + c(2, -2, -2, 2) / 7 * sqrt(6 / 5)))
  w <- switch(n, 2, c(1, 1), c(5, 8, 5) / 9,
    (18 + c(-1, 1, 1, -1) * sqrt(30)) / 36)
  o <- order(x)
  list(nodes = (x[o] + 1) / 2, weights = w[o] / 2)
}

## Nodes and log weights for the integral over u, per observation. The range
## (0, Inf) is cut at the u where either marginal's log-odds cross a level in
## .COP_LEV; each piece is weighed by the integrand at its ends, dropped if its
## mass is negligible, and otherwise split into enough equal sub-pieces in
## log u that m-point Gauss-Legendre resolves it. Returns long vectors with
## the observation index `row`, and an offset per row for the log-sum.
.cop_quad <- function(eps, S, su, sv, av, udist, vdist, copula, cpar,
                      tol = 1e-10, uA = .cop_u_at_lev(su, udist), m = 8L) {
  n <- length(eps)
  L <- .COP_LEV; nl <- length(L)
  vB <- .cop_v_at_lev(sv, av, vdist)
  UB <- outer(-eps, S * vB, "+")
  ## Where the noise at u = 0 is already past the last level, its density
  ## falls by about one log unit per unit of log-odds, so the same steps are
  ## taken from where it starts.
  B0 <- .cop_v_lo(S * eps, sv, av, vdist)
  far <- which(abs(B0) > 30)
  RL <- matrix(NA_real_, n, length(.COP_REL))
  if (length(far)) RL[far, ] <- outer(B0[far], S * .COP_REL, "+")
  UR <- matrix(Inf, n, length(.COP_REL))
  if (length(far)) UR[far, ] <- S * .cop_v_at_lev(sv, av, vdist, RL[far, ]) - eps[far]
  UB <- cbind(UB, UR)
  BR <- cbind(matrix(uA, n, nl, byrow = TRUE), UB)
  Aval <- cbind(matrix(L, n, nl, byrow = TRUE), matrix(NA_real_, n, ncol(UB)))
  Bval <- cbind(matrix(NA_real_, n, nl), matrix(L, n, nl, byrow = TRUE), RL)
  pos <- is.finite(UB) & UB > 0
  Aval[, nl + seq_len(ncol(UB))][pos] <- .cop_u_lo(UB[pos], su, udist)
  Bval[, seq_len(nl)] <- .cop_v_lo(S * (eps + BR[, seq_len(nl)]), sv, av, vdist)
  BR[!(BR > 0)] <- Inf
  K <- ncol(BR)
  rr <- rep(seq_len(n), K)
  o <- order(rr, as.vector(BR))
  BR <- matrix(as.vector(BR)[o], n, K, byrow = TRUE)
  Aval <- pmin(pmax(matrix(as.vector(Aval)[o], n, K, byrow = TRUE), -1e100), 1e100)
  Bval <- pmin(pmax(matrix(as.vector(Bval)[o], n, K, byrow = TRUE), -1e100), 1e100)
  fin <- is.finite(BR)
  ## Scout: the integrand, per unit of log u, at every breakpoint.
  sc <- matrix(-Inf, n, K)
  ld <- .cop_node_ld(BR[fin], matrix(eps, n, K)[fin], S, su, sv, av, udist,
    vdist, copula, cpar, parts = TRUE)
  if (is.null(ld)) return(NULL)
  sc[fin] <- ld$lv + ld$lu + ld$lc + log(BR[fin])
  sc[is.nan(sc)] <- -Inf
  ## Each component at each breakpoint, for the variation inside a piece.
  cmp <- lapply(list(ld$lv, ld$lu, ld$lc, log(BR[fin])), function(x) {
    out <- matrix(0, n, K); out[fin] <- x; out[!is.finite(out)] <- 0; out })
  ## Interior pieces.
  j1 <- seq_len(K - 1L); j2 <- j1 + 1L
  ok <- fin[, j2, drop = FALSE]
  lo <- log(BR[, j1, drop = FALSE]); wd <- log(BR[, j2, drop = FALSE]) - lo
  ok <- ok & wd > 0
  A1 <- Aval[, j1, drop = FALSE]; A2 <- Aval[, j2, drop = FALSE]
  B1 <- Bval[, j1, drop = FALSE]; B2 <- Bval[, j2, drop = FALSE]
  dAB <- abs(A2 - A1) + abs(B2 - B1)
  ## A copula ridge lies along A = B or A = -B; a piece that crosses one can
  ## hide a peak that neither end shows.
  cross <- (A1 - B1) * (A2 - B2) <= 0 | (A1 + B1) * (A2 + B2) <= 0 |
    pmin(abs(A1 - B1), abs(A2 - B2), abs(A1 + B1), abs(A2 + B2)) < 2
  ## Monotone components that move in opposite directions can peak inside a
  ## piece by up to the smaller of their total rise and total fall.
  up <- 0; dn <- 0
  for (x in cmp) {
    d <- x[, j2, drop = FALSE] - x[, j1, drop = FALSE]
    up <- up + pmax(d, 0); dn <- dn - pmin(d, 0)
  }
  lam <- abs(sc[, j2, drop = FALSE] - sc[, j1, drop = FALSE]) +
    2 * pmin(up, dn) + dAB * .cop_kappa(copula, cpar) * cross
  lam[!is.finite(lam)] <- 80
  score <- pmax(sc[, j1, drop = FALSE], sc[, j2, drop = FALSE]) + log(wd)
  score[!ok] <- -Inf
  ## Head piece (0, b_1) and tail piece (b_last, Inf); both in the deep tail
  ## of F_U, so a single piece each.
  nfin <- rowSums(fin)
  b1 <- BR[, 1L]
  bL <- BR[cbind(seq_len(n), nfin)]
  sL <- sc[cbind(seq_len(n), nfin)]
  hL <- pmin(if (identical(udist, "exponential")) rep_len(su, n) else su^2 / bL, sv)
  M <- pmax(apply(score, 1L, max), sc[, 1L], sL + log(hL / bL))
  M[!is.finite(M)] <- 0
  ## Prune, then cut. m-point Gauss-Legendre on exp(c x), x in (-1, 1), has
  ## relative error 2^(2m+1) (m!)^4 c^(2m+1) / ((2m+1) ((2m)!)^3). Each piece
  ## gets the rule and number of equal sub-pieces that bring that below tol
  ## over its share of the row's mass at the fewest nodes, so a piece holding
  ## 1e-20 of the mass costs two nodes.
  rowk <- (seq_along(score) - 1L) %% n + 1L
  lp <- pmin(score - M[rowk], 0)
  keep <- which(ok & lp + lam - abs(sc[, j2, drop = FALSE] - sc[, j1, drop = FALSE]) > -40)
  lp <- lp[keep]; lmk <- lam[keep]
  ms <- sort(unique(c(2L, 3L, 5L, m)))
  cost <- vapply(ms, function(mm) {
    lk <- log(2 * mm + 1) + 3 * lfactorial(2 * mm) - (2 * mm + 1) * log(2) - 4 * lfactorial(mm)
    Lam <- 2 * exp((log(tol) - lp + lk) / (2 * mm + 1))
    pmin(pmax(ceiling(lmk / Lam), 1), 1e6) * mm
  }, numeric(length(keep)))
  cost <- matrix(cost, ncol = length(ms))
  pick <- max.col(-cost, ties.method = "first")
  kk <- pmin(cost[cbind(seq_along(keep), pick)] / ms[pick], 64)
  row <- rowk[keep]
  Ul <- list(); LWl <- list(); Rl <- list()
  for (g in unique(pick)) {
    sel <- which(pick == g)
    mm <- ms[g]
    gl <- .cop_gl01(mm)
    sub <- rep(sel, kk[sel])
    hw <- wd[keep][sub] / kk[sub]
    slo <- lo[keep][sub] + (sequence(kk[sel]) - 1) * hw
    T <- outer(hw, gl$nodes) + slo
    Ul[[g]] <- as.vector(exp(T))
    LWl[[g]] <- as.vector(outer(log(hw), log(gl$weights), "+") + T)
    Rl[[g]] <- rep(row[sub], mm)
  }
  U <- unlist(Ul); LW <- unlist(LWl); R <- unlist(Rl)
  gl <- .gauss_legendre_01(as.integer(m))
  ## Head and tail pieces: m nodes each.
  Uh <- outer(b1, gl$nodes); LWh <- log(outer(b1, gl$weights))
  tt <- gl$nodes / (1 - gl$nodes)
  Ut <- bL + outer(hL, tt); LWt <- log(outer(hL, gl$weights / (1 - gl$nodes)^2))
  list(U = c(as.vector(U), as.vector(Uh), as.vector(Ut)),
    LW = c(as.vector(LW), as.vector(LWh), as.vector(LWt)),
    row = c(R, rep(seq_len(n), 2L * m)), M = M)
}

## log f_eps by .cop_quad(); with `parts = TRUE`, also the nodes and their
## normalized weights, for the efficiency predictors.
.cop_log_dens <- function(eps, S, su, sv, av, udist, vdist, copula, cpar,
                          tol = 1e-10, uA = .cop_u_at_lev(su, udist), parts = FALSE) {
  n <- length(eps)
  q <- .cop_quad(eps, S, su, sv, av, udist, vdist, copula, cpar, tol, uA)
  if (is.null(q)) return(NULL)
  ld <- .cop_node_ld(q$U, eps[q$row], S, su, sv, av, udist, vdist, copula, cpar)
  if (is.null(ld)) return(NULL)
  lg <- ld + q$LW
  lg[is.nan(lg)] <- -Inf
  s <- rowsum(exp(lg - q$M[q$row]), q$row, reorder = TRUE)[, 1L]
  out <- log(s) + q$M
  if (!parts) return(out)
  list(logf = out, U = q$U, row = q$row, w = exp(lg - out[q$row]))
}

copsfm <- function(formula,
                   data,
                   copula = c("gaussian", "fgm", "frank",
                              "clayton", "clayton90", "clayton180", "clayton270",
                              "gumbel", "gumbel90", "gumbel180", "gumbel270",
                              "joe", "joe90", "joe180", "joe270",
                              "independent"),
                   inefdec = TRUE,
                   udist = c("hnormal", "exponential"),
                   vdist = c("normal", "logistic", "glogistic", "laplace"),
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
  .check_optim_method(Method)
  call <- match.call()
  data_supplied <- if (!missing(data)) data
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

  ## The integral over u is by Gauss-Legendre on pieces placed per
  ## observation (.cop_quad(), A74). n_nodes sets the target error in log f:
  ## 1e-10 at the default 128, tenfold looser or tighter per 32 nodes.
  q_tol <- min(max(1e-10 * 10^((128 - n_nodes) / 32), 1e-14), 1e-6)
  cA <- .cop_u_at_lev(1, udist)
  .log_dens <- function(eps, su, sv, av, cpar, parts = FALSE) {
    out <- .cop_log_dens(eps, S, su, sv, av, udist, vdist,
      if (has_cop) copula else "independent", cpar, q_tol, uA = su * cA,
      parts = parts)
    if (is.null(out)) return(if (parts) NULL else rep(NA_real_, length(eps)))
    out
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

  if (optHessian == FALSE & PSopt == FALSE) { opt <- .as_optim(bob1); st_err <- rep(NA, length(opt$par)) }
  if (optHessian == FALSE & PSopt == TRUE)  { opt <- opt00; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == TRUE) {
    ## The colMeans test catches an all-zero Hessian only; a singular but
    ## non-zero one made solve() error and discard a finished fit -- minutes of
    ## one, this being simulated ML (A65).
    st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA, length(opt$par))
    } else {
      tryCatch(suppressWarnings(sqrt(diag(solve(opt$hessian)))),
        error = function(e) rep(NA_real_, length(opt$par))
      )
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

  ## Reported on the reported scale and names: the delta-method factors above
  ## carry a missing estimation-scale standard error straight through (A65).
  .sm <- if (optHessian == TRUE) {
    .se_missing_report(out["st_err", ], colnames(out), opt$hessian,
      paste0("copsfm(copula = \"", copula, "\"): ")
    )
  } else {
    NULL
  }
  if (!is.null(.sm)) warning(.sm$message, call. = FALSE)

  ## The efficiency predictors, by the same quadrature -- and from the same
  ## closure -- that the likelihood used.
  eps <- S * as.numeric(y - X %*% th[i_b])
  qp <- .log_dens(eps, su, sv, av, if (has_cop) cpar else NULL, parts = TRUE)
  jlms <- as.numeric(rowsum(qp$w * qp$U, qp$row, reorder = TRUE))
  ## Battese-Coelli, E[exp(-u) | eps]. Free here, and it is the quantity the
  ## copula literature reports; exp(-jlms) is kept as `efficiency` because that
  ## is what this function has always returned under that name.
  exp_u_hat <- as.numeric(rowsum(qp$w * exp(-qp$U), qp$row, reorder = TRUE))

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
  results$inefdec <- isTRUE(inefdec)
  ## Appended, not rebuilt: re-listing twenty-one values against twenty-one
  ## names to add one field is how a misalignment gets introduced.
  if (!is.null(.sm)) results$se_unresolved <- .sm
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
  .sfa_keep_data(results, data_supplied)
}
