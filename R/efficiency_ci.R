## Horrace and Schmidt (1996) intervals for individual inefficiency


## Conditional on the fitted parameters, u_i | e_i is normal with mean mu_star
## and standard deviation sigma_star, truncated below at zero.
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

  ## Taken from the fit where it carries them, so this and print()/summary()
  ## cannot drift. NHN and NHN_Z never store u_hat, so that one is usually
  ## rebuilt -- from the same posterior, so it is the same estimator.
  u_hat <- object$u_hat
  if (is.null(u_hat)) {
    u_hat <- .jlms_u(post$mu_star, post$sigma_star)
  }
  te_hat <- object$exp_u_hat
  n <- length(ci$lower)
  ## Each stored predictor is a fixed function of the posterior, so recompute
  ## it and compare. Length alone does not settle it: a vector of the right
  ## length and the wrong order passes a length check and is still wrong (A69).
  expected <- list(
    u_hat = function() as.numeric(.jlms_u(post$mu_star, post$sigma_star)),
    exp_u_hat = function() {
      as.numeric(.te_battese_coelli(post$mu_star, post$sigma_star))
    }
  )
  for (nm in c("u_hat", "exp_u_hat")) {
    v <- object[[nm]]
    if (is.null(v)) next
    ## A data frame column recycles a vector whose length divides the row
    ## count, so a mismatched predictor is misaligned silently (issue #2).
    if (length(v) != n) {
      stop("efficiency_ci(): `object$", nm, "` has ", length(v), " values but ",
        "the fitted posterior has ", n, ", so they cannot be matched to the ",
        "same observations.",
        call. = FALSE
      )
    }
    want <- expected[[nm]]()
    gap <- max(abs(as.numeric(v) - want))
    if (!is.finite(gap) || gap > 1e-8 * max(1, max(abs(want)))) {
      stop("efficiency_ci(): `object$", nm, "` is not what this fit's ",
        "posterior gives -- they differ by ", format(gap, digits = 3),
        ", so the stored values do not belong to these observations. Refit, ",
        "or remove the field to have it rebuilt from the posterior.",
        call. = FALSE
      )
    }
  }

  out <- data.frame(row.names = seq_along(ci$lower))
  if (type %in% c("both", "u")) {
    out$u_lower <- ci$lower
    out$u_hat <- as.numeric(u_hat)
    out$u_upper <- ci$upper
  }
  if (type %in% c("both", "te")) {
    ## Monotone decreasing transform, so the bounds swap.
    out$te_lower <- exp(-ci$upper)
    ## A fit stripped of exp_u_hat still carries the posterior that field was
    ## built from, so rebuild the same Battese-Coelli score rather than
    ## reporting exp(-JLMS), which is a different estimator (A68).
    out$te_hat <- if (is.null(te_hat)) {
      ## as.numeric(): mu_star arrives as a one-column matrix from some
      ## models, and a matrix column would change this field's type.
      as.numeric(.te_battese_coelli(post$mu_star, post$sigma_star))
    } else {
      as.numeric(te_hat)
    }
    out$te_upper <- exp(-ci$lower)
  }

  attr(out, "level") <- level
  out
}
