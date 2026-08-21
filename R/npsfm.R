## ---------------------------------------------------------------------------
## npsfm(): nonparametric stochastic frontier models
##
## Ported from Christopher Parmeter's research scripts (base code/Parmeter
## codes/npsfm*.R). Those scripts were drafts: the dispatcher passed
## unevaluated `match.call()` components where values were expected, several
## branches referenced undefined objects, and the helper routines relied on
## np's unexported internals. What is kept from them is the algebra and the
## paper references; the plumbing is written against this package's own
## conventions.
##
## Two estimators are implemented:
##
##   "FLW"  Fan, Li and Weersink (1996, JBES). Two steps. Fit E[y|x]
##          nonparametrically, then recover the scale parameters from the
##          residuals -- by maximizing Fan-Li-Weersink's concentrated
##          likelihood in lambda for the half-normal case, or by inverting
##          central moments for the exponential, gamma and uniform cases.
##          sigma_u and sigma_v are constants; only the frontier is smooth.
##
##   "SVKZ" Simar, Van Keilegom and Zelenyuk (2017, JPA). Local method of
##          moments. Three local-linear regressions -- y on x, then the
##          squared and cubed residuals on x -- give sigma_u(x) and
##          sigma_v(x) pointwise, so both variance components vary with the
##          covariates. No optimizer runs.
##
## Both need kernel regression from the np package, which is in Suggests
## rather than Imports: it is required only by this function, and 1.1.3
## deliberately reduced the package's hard dependencies.
## ---------------------------------------------------------------------------

.npsfm_need_np <- function() {
  if (!requireNamespace("np", quietly = TRUE)) {
    stop("npsfm() requires the 'np' package for kernel regression and ",
      "bandwidth selection, which is not installed. ",
      'Install it with install.packages("np").',
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## Fan, Li and Weersink's concentrated log-likelihood, their equation (14) and
## the discussion beneath (15). sigma is profiled out of the problem as a
## function of lambda, so only lambda is searched over.
##
## The original script maximized this with bbmle::mle2(). A one-dimensional
## bounded search does the same job without the extra dependency.
.flw_neg_cll <- function(lambda, e) {
  if (!is.finite(lambda) || lambda <= 0) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  sc <- sqrt(2 * lambda^2 / (pi * (1 + lambda^2)))
  s2 <- mean(e^2) / (1 - sc^2)
  if (!is.finite(s2) || s2 <= 0) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  si <- sqrt(s2)
  ## Mean-correct the residuals: E[eps] = -E[u] = -sigma*sc
  ep <- e - si * sc
  ll <- -length(ep) * log(si) +
    sum(pnorm(-ep * lambda / si, log.p = TRUE)) -
    0.5 * sum(ep^2) / s2
  if (!is.finite(ll)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -ll
}

## Kernel regression of ydat on xdat, returning fitted values, gradients and
## residuals. Deliberately uses np's xdat/ydat interface rather than its
## formula interface: the formula interface resolves `data` by non-standard
## evaluation in the caller's frame, which is what made the original
## FLW.fun() fail with "'data' must be a data.frame" when called from inside
## another function.
.npsfm_kreg <- function(xdat, ydat, regtype, bw.sel, bw = NULL, nmulti = 1L) {
  if (is.null(bw)) {
    bws <- np::npregbw(
      xdat = xdat, ydat = ydat, regtype = regtype,
      bwmethod = bw.sel, nmulti = nmulti
    )
  } else {
    bws <- np::npregbw(
      xdat = xdat, ydat = ydat, regtype = regtype,
      bandwidth.compute = FALSE
    )
    bws$bw <- bw
  }
  fit <- np::npreg(bws = bws, gradients = TRUE, residuals = TRUE)
  ## Namespace-qualify everything: np is only in Suggests, so it is not
  ## attached when npsfm() runs. fitted()/residuals() are stats generics that
  ## np provides methods for, but `gradients()` is np's OWN generic and is
  ## invisible unless np is on the search path.
  list(
    bws = bws,
    fitted = as.numeric(stats::fitted(fit)),
    grad = as.matrix(np::gradients(fit)),
    resid = as.numeric(stats::residuals(fit))
  )
}


## ---------------------------------------------------------------------------
## FLW: Fan, Li and Weersink (1996)
## ---------------------------------------------------------------------------
.npsfm_flw <- function(xdat, ydat, dist, regtype, bw.sel, bw, cost) {
  cost_sc <- if (cost) -1 else 1
  k <- .npsfm_kreg(
    xdat = xdat, ydat = cost_sc * ydat, regtype = regtype,
    bw.sel = bw.sel, bw = bw
  )
  e <- k$resid
  m2 <- mean(e^2)
  m3 <- mean(e^3)
  m4 <- mean(e^4)

  sigma_u <- sigma_v <- lambda <- sigma <- NA_real_
  theta <- bhat <- NA_real_

  if (dist == "hn") {
    ## Identification of sigma_u rests on the residuals being negatively
    ## skewed. If they are not, the concentrated likelihood is maximized at
    ## lambda -> 0 and sigma_u collapses to essentially zero. That is the
    ## correct answer to the question asked, but it is silent, and the usual
    ## cause is a misspecified orientation -- fitting production data with
    ## cost = TRUE, or the reverse -- rather than a genuinely efficient
    ## sample. Say so.
    if (m3 >= 0) {
      warning("FLW: residuals are not negatively skewed (third moment ",
        signif(m3, 3), "), so sigma_u is not identified and will be ",
        "driven to zero. Check the `cost` orientation.",
        call. = FALSE
      )
    }
    ## Search lambda on a wide bounded interval; the concentrated likelihood
    ## is one-dimensional so this is both cheap and reliable.
    opt <- stats::optimize(.flw_neg_cll, interval = c(1e-4, 50), e = e)
    lambda <- opt$minimum
    sc <- sqrt(2 * lambda^2 / (pi * (1 + lambda^2)))
    sigma <- sqrt(m2 / (1 - sc^2))
    sigma_v <- sigma / sqrt(1 + lambda^2)
    sigma_u <- sigma_v * lambda
    mu <- sigma * sc
  } else if (dist == "exp") {
    ## u ~ Exp(rate theta): mu3(eps) = -2/theta^3
    if (m3 >= 0) {
      warning("FLW: third residual moment is non-negative (wrong skew); ",
        "the exponential scale is not identified from these data.",
        call. = FALSE
      )
      theta <- NA_real_
      mu <- 0
    } else {
      theta <- (-2 / m3)^(1 / 3)
      sigma_u <- 1 / theta
      sigma_v <- sqrt(max(m2 - 1 / theta^2, 0))
      mu <- 1 / theta
    }
  } else if (dist == "gamma") {
    ## u ~ Gamma(P, theta); uses the third and fourth central moments
    m4star <- m4 - 3 * m2^2
    theta <- -3 * m3 / m4star
    p_hat <- -theta^3 * m3 / 2
    sigma_v <- sqrt(max(m2 - p_hat / theta^2, 0))
    sigma_u <- if (is.finite(p_hat) && p_hat > 0) sqrt(p_hat) / theta else NA_real_
    mu <- p_hat / theta
  } else if (dist == "unif") {
    ## u ~ U(0, b)
    bhat2 <- (3 * m2^2 - m4) * 120
    if (!is.finite(bhat2) || bhat2 < 0) {
      warning("FLW: implied uniform upper bound is negative; setting it to zero.",
        call. = FALSE
      )
      bhat <- 0
    } else {
      bhat <- bhat2^(1 / 4)
    }
    sigma_v <- sqrt(max(m2 - bhat^2 / 12, 0))
    sigma_u <- bhat / sqrt(12)
    mu <- bhat / 2
  }

  ## The kernel fit estimates E[y|x] = m(x) - E[u]; shift it back up to the
  ## frontier. Gradients are unaffected because the shift is constant.
  frontier <- cost_sc * (k$fitted + mu)
  frontier_grad <- cost_sc * k$grad
  eps <- cost_sc * e

  list(
    frontier = frontier, frontier_grad = frontier_grad, residuals = eps,
    sigma_u = sigma_u, sigma_v = sigma_v, lambda = lambda, sigma = sigma,
    theta = theta, b = bhat, mean_correction = mu,
    conditional_mean = cost_sc * k$fitted, bws = k$bws
  )
}


## ---------------------------------------------------------------------------
## SVKZ: Simar, Van Keilegom and Zelenyuk (2017)
## ---------------------------------------------------------------------------
.npsfm_svkz <- function(xdat, ydat, bw.sel, cost) {
  cost_sc <- if (cost) -1 else 1
  yy <- cost_sc * ydat

  ## r1: conditional mean
  r1 <- .npsfm_kreg(xdat = xdat, ydat = yy, regtype = "ll", bw.sel = bw.sel)
  e1 <- r1$resid

  ## r2, r3: conditional second and third moments of the residual
  r2 <- .npsfm_kreg(xdat = xdat, ydat = e1^2, regtype = "ll", bw.sel = bw.sel)
  r3 <- .npsfm_kreg(xdat = xdat, ydat = e1^3, regtype = "ll", bw.sel = bw.sel)

  ## SVKZ (3.5): sigma_u(x)^3 = sqrt(pi/2) * (pi/(pi-4)) * r3(x).
  ## pi/(pi-4) is negative, so a correctly signed (negative) r3 gives a
  ## positive cube. Wrong-skew points are floored at zero by (3.7).
  sigu3 <- sqrt(pi / 2) * (pi / (pi - 4)) * r3$fitted
  wrong_skew <- sigu3 <= 0
  sigma_u <- pmax(sigu3, 0)^(1 / 3)

  ## (3.8): sigma_v(x)^2 = r2(x) - sigma_u(x)^2 (pi-2)/pi, floored at zero
  sigma_v2 <- pmax(r2$fitted - sigma_u^2 * (pi - 2) / pi, 0)

  ## (3.2)/(3.9): mean correction E[u|x] = sqrt(2/pi) sigma_u(x)
  mu <- sqrt(2 / pi) * sigma_u

  ## Gradient of the correction, by the chain rule through (3.5):
  ##   d sigma_u/dx = (1/3) * sigu3^(-2/3) * d(sigu3)/dx
  ## The source script wrote b3hat^(-2/3) here, which differs from
  ## sigu3^(-2/3) * b3hat except in the degenerate case, and is dimensionally
  ## wrong once x has more than one column.
  dsigu3 <- sqrt(pi / 2) * (pi / (pi - 4)) * r3$grad
  scale_g <- (1 / 3) * ifelse(wrong_skew, 0, pmax(sigu3, .SFA_CONSTANTS$MIN_POSITIVE)^(-2 / 3))
  sigma_u_grad <- dsigu3 * scale_g
  sigma_u_grad[wrong_skew, ] <- 0

  frontier <- cost_sc * (r1$fitted + mu)
  frontier_grad <- cost_sc * (r1$grad + sqrt(2 / pi) * sigma_u_grad)

  list(
    frontier = frontier, frontier_grad = frontier_grad,
    residuals = cost_sc * e1,
    sigma_u = sigma_u, sigma_v = sqrt(sigma_v2),
    sigma_u_grad = sigma_u_grad, mean_correction = mu,
    conditional_mean = cost_sc * r1$fitted, wrong_skew = wrong_skew,
    bws = list(r1 = r1$bws, r2 = r2$bws, r3 = r3$bws)
  )
}


## ---------------------------------------------------------------------------
## User-facing entry point
## ---------------------------------------------------------------------------
npsfm <- function(formula, data, method = c("FLW", "SVKZ"),
                  dist = c("hn", "exp", "gamma", "unif"),
                  regtype = c("lc", "ll"), bw.sel = c("cv.ls", "cv.aic"),
                  bw = NULL, cost = FALSE, eff = TRUE) {
  .npsfm_need_np()

  method <- match.arg(method)
  regtype <- match.arg(regtype)
  bw.sel <- match.arg(bw.sel)
  dist <- match.arg(dist)

  if (method == "SVKZ" && dist != "hn") {
    stop("method = \"SVKZ\" is derived for the normal-half normal model only; ",
      "dist = \"", dist, "\" is not available for it.",
      call. = FALSE
    )
  }

  cl <- match.call()

  ## The frontier is the first pipe segment; npsfm() has no variance
  ## determinants of its own, so any further segments are an error rather
  ## than something to be silently dropped.
  parsed <- .parse_pipe_formula(formula)
  if (!is.null(parsed$formula_z)) {
    stop("npsfm() takes a single-part formula (y ~ x1 + x2). ",
      "Variance determinants enter nonparametrically through the ",
      "covariates themselves; there is no `| z` segment.",
      call. = FALSE
    )
  }

  mf <- stats::model.frame(parsed$formula_x, data = data, na.action = stats::na.omit)
  ydat <- as.numeric(stats::model.response(mf))
  xdat <- mf[, -1, drop = FALSE]
  if (!ncol(xdat)) {
    stop("npsfm() needs at least one covariate on the right-hand side.", call. = FALSE)
  }
  n_obs <- length(ydat)

  start_time <- Sys.time()
  fit <- if (method == "FLW") {
    .npsfm_flw(xdat, ydat, dist, regtype, bw.sel, bw, cost)
  } else {
    .npsfm_svkz(xdat, ydat, bw.sel, cost)
  }
  total_time <- Sys.time() - start_time

  ## Observation-level inefficiency predictions. The composed residual is
  ## measured against the corrected frontier, so eps = y - frontier.
  eps <- ydat - fit$frontier
  u_hat <- exp_u_hat <- NULL
  if (isTRUE(eff)) {
    sc <- if (cost) -1 else 1
    su <- fit$sigma_u
    sv <- fit$sigma_v
    if (dist == "hn" && all(is.finite(c(su, sv))) && all(sv > 0)) {
      s2 <- su^2 + sv^2
      mu_star <- -sc * eps * su^2 / s2
      sig_star <- su * sv / sqrt(s2)
      u_hat <- .jlms_u(mu_star, sig_star)
      exp_u_hat <- .te_battese_coelli(mu_star, sig_star)
    } else if (dist == "exp" && is.finite(fit$theta) && is.finite(sv) && sv > 0) {
      mu_star <- -sc * eps - sv^2 * fit$theta
      u_hat <- .jlms_u(mu_star, rep(sv, length(mu_star)))
      exp_u_hat <- .te_battese_coelli(mu_star, rep(sv, length(mu_star)))
    }
  }

  out <- list(
    method = method, dist = dist, formula = formula, call = cl,
    frontier = fit$frontier, frontier.grad = fit$frontier_grad,
    conditional.mean = fit$conditional_mean,
    residuals = eps, mean.correction = fit$mean_correction,
    sigma.u = fit$sigma_u, sigma.v = fit$sigma_v,
    lambda = fit$lambda, sigma = fit$sigma,
    theta = fit$theta, b = fit$b,
    sigma.u.grad = fit$sigma_u_grad, wrong.skew = fit$wrong_skew,
    u_hat = u_hat, exp_u_hat = exp_u_hat,
    bws = fit$bws, cost = cost, regtype = regtype, bw.sel = bw.sel,
    nobs = n_obs, total_time = total_time, data = data
  )
  out <- out[!vapply(out, is.null, logical(1))]
  class(out) <- "npsfareg"
  out
}


## ---------------------------------------------------------------------------
## Methods. "npsfareg" is deliberately a separate class from "sfareg": there
## is no parameter vector with standard errors here, so coef()/vcov()/logLik()
## have nothing to return and it would be misleading to inherit them.
## ---------------------------------------------------------------------------
print.npsfareg <- function(x, ...) {
  cat("Nonparametric stochastic frontier model\n")
  cat("Method:       ", x$method,
    switch(x$method,
      FLW = " (Fan, Li and Weersink 1996)",
      SVKZ = " (Simar, Van Keilegom and Zelenyuk 2017)"
    ), "\n",
    sep = ""
  )
  cat("Distribution: ", x$dist, "\n", sep = "")
  cat("Frontier:     ", if (x$cost) "cost" else "production", "\n", sep = "")
  cat("Observations: ", x$nobs, "\n", sep = "")
  if (length(x$sigma.u) == 1L) {
    cat("\nsigma_u = ", signif(x$sigma.u, 5),
      "   sigma_v = ", signif(x$sigma.v, 5), "\n",
      sep = ""
    )
    if (!is.null(x$lambda) && is.finite(x$lambda)) {
      cat("lambda  = ", signif(x$lambda, 5),
        "   sigma   = ", signif(x$sigma, 5), "\n",
        sep = ""
      )
    }
  } else {
    cat("\nsigma_u(x) and sigma_v(x) vary by observation:\n")
    print(summary(data.frame(sigma_u = x$sigma.u, sigma_v = x$sigma.v)))
    if (!is.null(x$wrong.skew) && any(x$wrong.skew)) {
      cat("\n", sum(x$wrong.skew), " of ", x$nobs,
        " points had the wrong residual skew and were floored at sigma_u = 0.\n",
        sep = ""
      )
    }
  }
  if (!is.null(x$exp_u_hat)) {
    cat("\nMean predicted efficiency E[exp(-u)|e]: ",
      signif(mean(x$exp_u_hat, na.rm = TRUE), 4), "\n",
      sep = ""
    )
  }
  invisible(x)
}

summary.npsfareg <- function(object, ...) {
  print(object)
  cat("\nCall:\n")
  print(object$call)
  cat("\nFrontier fit:\n")
  print(summary(object$frontier))
  cat("\nComposed residuals (y - frontier):\n")
  print(summary(object$residuals))
  invisible(object)
}

fitted.npsfareg <- function(object, ...) object$frontier

residuals.npsfareg <- function(object, ...) object$residuals

nobs.npsfareg <- function(object, ...) object$nobs
