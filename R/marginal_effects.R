## ---------------------------------------------------------------------------
## Marginal effects of the variance determinants on inefficiency
##
## The _Z models let sigma_u depend on covariates, and sfa reported the
## delta coefficients and nothing else. A delta is not what applied work
## reports, because it is not interpretable on its own: it is a coefficient in
## a log-link for a scale parameter, so its units are neither those of u nor
## those of z, and its magnitude cannot be compared across models whose link
## differs. What gets reported is d E[u]/d z_k and d Var[u]/d z_k.
##
## Both are closed form. Writing s = sigma_u(z) and taking the half-normal
## case, E[u] = s sqrt(2/pi) and Var[u] = s^2 (1 - 2/pi), so every derivative
## reduces to d s/d z_k times a constant, and the constant cancels into E[u]
## or Var[u] itself:
##
##   SD link,        sigma_u = exp(z'delta):        d s/d z_k = delta_k s
##   variance link,  sigma_u = sqrt(exp(z'delta)):  d s/d z_k = (delta_k/2) s
##
## giving, for BOTH the half-normal and the exponential,
##
##   SD link:        dE[u]/dz_k = delta_k E[u],       dVar[u]/dz_k = 2 delta_k Var[u]
##   variance link:  dE[u]/dz_k = (delta_k/2) E[u],   dVar[u]/dz_k = delta_k Var[u]
##
## THE FACTOR OF TWO IS THE POINT. sfm()'s NHN_Z/NE_Z put the linear predictor
## on the standard deviation and psfm()'s TRE_Z/GTRE_Z put it on the variance
## (see CLAUDE.md, and entry C1 of the horserace gap list). Reading a delta
## from one family with the other family's convention in mind is off by
## exactly this factor. Reporting the marginal effect rather than the
## coefficient removes the trap, because the marginal effect is on the scale
## of u either way.
##
## The effects are per observation, because sigma_u varies with z. The
## average over observations -- the "average marginal effect" -- is the
## conventional summary and is returned as an attribute.
## ---------------------------------------------------------------------------


## marginal_effects(object, average = FALSE)
##
## Constant columns of the z design are dropped: the derivative with respect to
## an intercept is not a marginal effect, and reporting delta_0 * E[u] in a
## column headed "(Intercept)" invites it to be read as one.
##
## Standard errors are NOT reported. The marginal effect is a nonlinear
## function of both delta and z, so its sampling distribution needs either the
## delta method through the full vcov or a bootstrap; neither is free, and a
## number presented without one would be read as inference. psfm_bootstrap()
## is the existing route for that on the panel side.
marginal_effects <- function(object, average = FALSE) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit, as returned by sfm().", call. = FALSE)
  }
  zs <- object$z_spec
  if (is.null(zs)) {
    stop(
      "marginal_effects() needs a model whose inefficiency scale depends on ",
      "covariates. This fit is model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\", which has a single homoskedastic sigma_u and therefore no ",
      "marginal effect to report. Use one of the _Z models with a formula ",
      "such as y ~ x1 + x2 | z.",
      call. = FALSE
    )
  }

  Z <- zs$Z
  delta <- zs$delta
  eta <- as.numeric(Z %*% delta)

  ## sigma_u on the scale each family actually uses.
  sigma_u <- if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))

  ## Moments of u given sigma_u.
  if (identical(zs$family, "halfnormal")) {
    e_u <- sigma_u * sqrt(2 / pi)
    var_u <- sigma_u^2 * (1 - 2 / pi)
  } else if (identical(zs$family, "exponential")) {
    e_u <- sigma_u
    var_u <- sigma_u^2
  } else {
    stop("unsupported inefficiency family: ", zs$family, call. = FALSE)
  }

  ## d sigma_u / d z_k = scale_k * sigma_u, with scale_k = delta_k under the
  ## SD link and delta_k/2 under the variance link.
  half <- if (identical(zs$link, "sd")) 1 else 0.5

  keep <- vapply(seq_len(ncol(Z)), function(j) length(unique(Z[, j])) > 1L, logical(1))
  if (!any(keep)) {
    stop("every column of the variance-determinant design is constant, so ",
      "there is no covariate to differentiate with respect to.",
      call. = FALSE
    )
  }
  nm <- colnames(Z)[keep]
  d <- delta[keep]

  me_e <- outer(e_u, half * d)
  me_v <- outer(var_u, 2 * half * d)
  colnames(me_e) <- paste0("dE_u.d", nm)
  colnames(me_v) <- paste0("dVar_u.d", nm)

  out <- data.frame(sigma_u = sigma_u, E_u = e_u, Var_u = var_u,
                    me_e, me_v, check.names = FALSE)

  avg <- colMeans(cbind(me_e, me_v))
  attr(out, "average") <- avg
  attr(out, "link") <- zs$link
  attr(out, "family") <- zs$family

  if (isTRUE(average)) {
    return(avg)
  }
  out
}
