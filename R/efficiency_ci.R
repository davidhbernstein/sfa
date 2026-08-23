## ---------------------------------------------------------------------------
## Horrace and Schmidt (1996) intervals for individual inefficiency
##
## sfa has always reported point predictions of u_i -- u_hat (Jondrow et al.
## 1982) and exp_u_hat (Battese and Coelli 1988) -- and nothing about how
## sharply either is pinned down. Both are posterior MEANS of u given the
## composed residual e_i, and the posterior they average over is available in
## closed form, so the corresponding interval costs no extra estimation. This
## is the standard thing applied work reports alongside an efficiency ranking,
## and a reviewer can reasonably ask why a ranking is presented without it.
##
## The arithmetic lives in .horrace_schmidt_ci() (matrix_utils.R), next to the
## two point predictors it shares (mu_star, sigma_star) with. This file is the
## user-facing wrapper: it finds those quantities on the fit, checks the model
## is one whose posterior is actually a truncated normal, and lays the result
## out as a data frame.
## ---------------------------------------------------------------------------


## efficiency_ci(object, level = 0.95, type = c("both", "u", "te"))
##
## Conditional on the fitted parameters, u_i | e_i is normal with mean mu_star
## and standard deviation sigma_star, truncated below at zero. Inverting that
## truncated normal at (1 - level)/2 and 1 - (1 - level)/2 gives the bounds in
## closed form; the efficiency bounds follow by monotonicity, with the
## endpoints swapped, since exp(-u) is decreasing in u.
##
## THESE INTERVALS CONDITION ON THE ESTIMATED PARAMETERS. They say where u_i
## sits given e_i and a KNOWN frontier, and make no allowance for sampling
## error in the slopes or the variance parameters -- which is why they do not
## narrow as n grows. They measure the irreducible difficulty of splitting one
## residual into noise and inefficiency, not estimation uncertainty. A wide
## interval is a warning against reading much into that unit's rank.
##
## Available for the models whose posterior really is a truncated normal:
## NHN, NHN_Z, NE and NTN. Every other model_name has a posterior of a
## different shape, for which these formulas do not hold, so this errors
## rather than returning a misleading number.
##
## References: Horrace and Schmidt (1996), Journal of Productivity Analysis 7,
## 257-282; Jondrow et al. (1982); Battese and Coelli (1988). Full citations
## in man/efficiency_ci.Rd.
efficiency_ci <- function(object, level = 0.95, type = c("both", "u", "te")) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit, as returned by sfm().", call. = FALSE)
  }
  type <- match.arg(type)

  post <- object$u_posterior
  if (is.null(post)) {
    stop(
      "efficiency_ci() needs the posterior of u given the composed residual, ",
      "which is a truncated normal only for model_name \"NHN\", \"NHN_Z\", ",
      "\"NE\" and \"NTN\". This fit is model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\", whose posterior has a different shape, so the Horrace-Schmidt ",
      "bounds do not apply to it.",
      call. = FALSE
    )
  }

  ci <- .horrace_schmidt_ci(post$mu_star, post$sigma_star, level = level)

  ## The point predictors are taken from the fit rather than recomputed, so
  ## what this function reports and what print()/summary() report cannot drift.
  u_hat <- object$u_hat
  if (is.null(u_hat)) {
    u_hat <- .jlms_u(post$mu_star, post$sigma_star)
  }
  te_hat <- object$exp_u_hat

  out <- data.frame(row.names = seq_along(ci$lower))
  if (type %in% c("both", "u")) {
    out$u_lower <- ci$lower
    out$u_hat <- as.numeric(u_hat)
    out$u_upper <- ci$upper
  }
  if (type %in% c("both", "te")) {
    ## Monotone decreasing transform, so the bounds swap.
    out$te_lower <- exp(-ci$upper)
    out$te_hat <- if (is.null(te_hat)) exp(-as.numeric(u_hat)) else as.numeric(te_hat)
    out$te_upper <- exp(-ci$lower)
  }

  attr(out, "level") <- level
  out
}
