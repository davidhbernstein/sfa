## Bonanno, De Giovanni and Domma (2017), "The 'wrong skewness' problem: a
## re-specification of stochastic frontiers", Journal of Productivity Analysis
## 47:49-64.  Gap L31.
##
## Their argument, and the reason this function exists: the composed error
## eps = v - u has a third central moment that splits into three additive
## pieces, and the classic model forces TWO of them to zero. Write
## Utilde = u - E u and Vtilde = v - E v; then
##
##   E[(eps - E eps)^3] = E[Vtilde^3] - E[Utilde^3] - 3E[Vtilde^2 Utilde]
##                                                  + 3E[Vtilde Utilde^2]
##
## and the last two vanish under independence. A normal v kills the first as
## well, leaving nothing but the inefficiency term -- which is why, under the
## textbook specification, a positively skewed residual has nowhere to go but
## sigma_u = 0. Skew in the noise, or dependence between the two components,
## are perfectly ordinary explanations that the textbook model cannot express.
## That is the sense in which the wrong-skewness anomaly is, in their phrase,
## an ill-posed problem.
##
## For an exponential u, a generalized-logistic v and an FGM copula the three
## pieces have the closed form of their Eq. (11); this function computes the
## dependence piece by quadrature instead, which agrees with that closed form
## (see tests/testthat/test-skewdecomp.R) and also covers the marginals and
## copulas Eq. (11) does not.

## Quantile functions on the (0,1) scale, so the cross moments become an
## integral over the unit square against the copula density.
.sd_qu <- function(a, su, udist) {
  if (identical(udist, "exponential")) -su * log1p(-a) else su * stats::qnorm((1 + a) / 2)
}
.sd_qv <- function(b, sv, av, vdist) {
  if (identical(vdist, "normal")) sv * stats::qnorm(b) else .gl_q(b, av, sv)
}

## Marginal central moments. The third-moment COEFFICIENTS already live in
## copsfm.R as .cop_u_k3() and .cop_v_v2() etc., and are reused so the two
## files cannot drift apart.
.sd_u_moments <- function(su, udist) {
  list(mean = .cop_u_m1(udist) * su, var = .cop_u_v2(udist) * su^2,
    m3 = -.cop_u_k3(udist) * su^3
  )
}
.sd_v_moments <- function(sv, av, vdist) {
  m3 <- if (identical(vdist, "normal")) 0 else {
    sv^3 * (psigamma(av, 2L) - psigamma(1, 2L))
  }
  ## NOT .cop_v_v2(), which is the ORDINARY logistic constant pi^2/3 and holds
  ## only at alpha_v = 1. The generalized logistic's variance moves with the
  ## shape: delta^2 [psi'(alpha) + psi'(1)].
  vr <- if (identical(vdist, "normal")) sv^2 else {
    sv^2 * (trigamma(av) + trigamma(1))
  }
  list(mean = 0, var = vr, m3 = m3)
}

## Gauss-Legendre on (0,1)^2. The copula density's FIRST argument is the noise
## margin and its second the inefficiency margin, matching copsfm(); the
## rotated families are not symmetric in their arguments, so the order is not
## a free choice. 160 nodes per side is where the FGM case stops moving in the
## ninth digit.
.sd_grid <- function(su, sv, av, udist, vdist, copula, cpar, nodes) {
  gl <- .gauss_legendre_01(as.integer(nodes))
  a <- gl$nodes
  w <- gl$weights
  Ut <- .sd_qu(a, su, udist) - .sd_u_moments(su, udist)$mean
  Vt <- .sd_qv(a, sv, av, vdist)
  if (identical(copula, "independent") || !is.finite(cpar)) {
    D <- outer(w, w)
  } else {
    C <- matrix(.cop_logc_rot(rep(a, times = nodes), rep(a, each = nodes),
      cpar, copula
    ), nrow = nodes)
    C[!is.finite(C)] <- -Inf
    D <- exp(C) * outer(w, w)
  }
  list(Ut = Ut, Vt = Vt, D = D, a = a, w = w)
}

.sd_cross <- function(g) {
  c(v2u = sum(outer(g$Vt^2, g$Ut) * g$D),
    vu2 = sum(outer(g$Vt, g$Ut^2) * g$D),
    vu  = sum(outer(g$Vt, g$Ut) * g$D)
  )
}

## The median of eps = v - u. No closed form exists for any of these pairs.
## P(eps <= q) = int_0^1 int_0^{G_V(q + Q_U(a))} c(b, a) db da: the inner
## integral is the copula density over the noise margin at fixed inefficiency,
## so no copula CDF is needed.
.sd_median <- function(su, sv, av, udist, vdist, copula, cpar, nodes) {
  gl <- .gauss_legendre_01(as.integer(nodes))
  a <- gl$nodes
  w <- gl$weights
  U <- .sd_qu(a, su, udist)
  hascop <- !identical(copula, "independent") && is.finite(cpar)
  gi <- .gauss_legendre_01(64L)
  Fe <- function(q) {
    B <- .cop_v_p(q + U, sv, av, vdist)
    B <- pmin(pmax(B, 0), 1)
    if (!hascop) return(sum(w * B))
    ## inner node b = B[j] * gi$nodes[i], weight B[j] * gi$weights[i]
    bb <- outer(gi$nodes, B)
    lc <- .cop_logc_rot(as.numeric(bb), rep(a, each = length(gi$nodes)),
      cpar, copula
    )
    M <- matrix(lc, nrow = length(gi$nodes))
    M[!is.finite(M)] <- -Inf
    inner <- colSums(exp(M) * gi$weights) * B
    sum(w * inner)
  }
  hi <- 40 * (su + sv) + 1
  stats::uniroot(function(q) Fe(q) - 0.5, c(-hi, hi), tol = 1e-8)$root
}

skewness_decomp <- function(object, nodes = 160L) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be a fitted stochastic frontier model of class ",
      "\"sfareg\". Got an object of class ",
      paste(class(object), collapse = "/"), ".",
      call. = FALSE
    )
  }
  p <- object$out[, "par"]
  mn <- object$model_name

  if (identical(mn, "COP")) {
    udist <- object$udist
    vdist <- object$vdist
    su <- unname(p[["sigma_u"]])
    sv <- unname(p[[if (identical(vdist, "normal")) "sigma_v" else "delta_v"]])
    av <- if (identical(vdist, "glogistic")) unname(object$alpha_v) else 1
    cop <- object$copula
    cpar <- object$copula_par
  } else if (mn %in% c("NHN", "NE")) {
    udist <- if (identical(mn, "NHN")) "hnormal" else "exponential"
    vdist <- "normal"
    ## NHN under maximum likelihood reports the (lambda, sigma)
    ## reparameterization; under "cols"/"acols"/"cmle", and NE always, the
    ## two scales are reported directly. Both have to be read.
    if (all(c("sigu", "sigv") %in% names(p))) {
      su <- unname(p[["sigu"]])
      sv <- unname(p[["sigv"]])
    } else if (all(c("lambda", "sigma") %in% names(p))) {
      lam <- unname(p[["lambda"]])
      sg <- unname(p[["sigma"]])
      sv <- sg / sqrt(1 + lam^2)
      su <- lam * sv
    } else {
      stop("skewness_decomp(): could not find the two scale parameters on ",
        "this fit. Expected \"sigu\"/\"sigv\" or \"lambda\"/\"sigma\" ",
        "among ", paste(names(p), collapse = ", "), ".",
        call. = FALSE
      )
    }
    av <- 1
    cop <- "independent"
    cpar <- NA_real_
  } else {
    stop("skewness_decomp() needs a model whose two error components have ",
      "separately identified third moments. That means a copsfm() fit, or an ",
      "sfm() fit with model_name \"NHN\" or \"NE\". Got model_name \"", mn, "\".",
      call. = FALSE
    )
  }

  um <- .sd_u_moments(su, udist)
  vm <- .sd_v_moments(sv, av, vdist)
  gg <- .sd_grid(su, sv, av, udist, vdist, cop, cpar, nodes)
  cr <- .sd_cross(gg)

  ## Under independence the two cross terms are exactly zero -- E[Vtilde] and
  ## E[Utilde] are both zero and the expectations factor -- so they are set,
  ## not integrated. Quadrature would put 1e-6 of its own noise on a quantity
  ## that is analytically nil, and readers would reasonably wonder what it was.
  indep <- identical(cop, "independent") || !is.finite(cpar)
  comp <- c(
    inefficiency = -um$m3,
    noise = vm$m3,
    dependence = if (indep) 0 else unname(3 * (cr[["vu2"]] - cr[["v2u"]]))
  )
  tot <- sum(comp)

  ## Zenga's median-based measure, their Eq. (13). It is the quantity whose
  ## SIGN says which way the density actually leans; the third central moment
  ## can disagree with it, and when it does the "wrong skewness" reading of a
  ## positive m3 is the wrong reading.
  e_mean <- vm$mean - um$mean
  e_var <- vm$var + um$var - if (indep) 0 else 2 * cr[["vu"]]
  e_med <- .sd_median(su, sv, av, udist, vdist, cop, cpar, nodes)
  med3 <- tot + (e_mean - e_med)^3 + 3 * (e_mean - e_med) * e_var

  res <- list(
    components = comp, total = tot, share = comp / tot,
    mean = e_mean, median = e_med, variance = e_var, median_based = med3,
    model_name = mn, udist = udist, vdist = vdist, copula = cop,
    copula_par = cpar, call = match.call()
  )
  class(res) <- "sfa_skew_decomp"
  res
}

print.sfa_skew_decomp <- function(x, ...) {
  cat("Third-moment decomposition of the composed error\n")
  cat("Bonanno, De Giovanni and Domma (2017), Eq. (11)\n\n")
  cat("  inefficiency: ", format(x$components[["inefficiency"]], digits = 5), "\n")
  cat("  noise:        ", format(x$components[["noise"]], digits = 5), "\n")
  cat("  dependence:   ", format(x$components[["dependence"]], digits = 5), "\n")
  cat("  ------------------------------\n")
  cat("  total m3:     ", format(x$total, digits = 5), "\n\n")
  cat("  median-based (Eq. 13): ", format(x$median_based, digits = 5), "\n")
  cat("  E[eps] ", format(x$mean, digits = 5),
    "  Me[eps] ", format(x$median, digits = 5),
    "  Var[eps] ", format(x$variance, digits = 5), "\n", sep = ""
  )
  if (sign(x$total) != sign(x$median_based)) {
    cat("\n  NOTE: the third central moment and the median-based measure ",
      "DISAGREE in\n  sign. The mean-based measure is the misleading one ",
      "here; see the paper.\n", sep = ""
    )
  }
  invisible(x)
}
