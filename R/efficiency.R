## efficiency(): the three point predictors on either scale of the dependent
## variable. See notes/code_history/efficiency.md.
##
## Two gaps closed at once, because they are the same decision seen twice: what
## exactly is being reported as "efficiency".
##
##  - H7. The package assumed throughout that the dependent variable is logged,
##    without saying so anywhere in a signature. `logDepVar` makes the
##    assumption visible and gives the level-scale alternative.
##  - G7. Only two of the three usual predictors were available. `type = "mode"`
##    adds the third.

## Mode of the posterior of u given the composed residual. For the models whose
## posterior is a truncated normal N+(mu_star, sigma_star^2) this is simply
## mu_star censored at zero: the untruncated normal's mode is its mean, and
## truncation to the positive half-line moves it to the boundary whenever the
## mean is negative.
##
## The mode is the ONLY one of the three predictors that can be exactly zero,
## and it is zero for every observation whose posterior mean is negative --
## which is most of the efficient tail. That is a property of the estimator,
## not a defect: the most likely value of u really is the boundary there.
.u_mode <- function(mu_star) pmax(mu_star, 0)

## The persistent efficiency of a four-component panel fit. psfm()'s GTRE,
## GTRE_FML and GTRE_Z store it as $H, one value per FIRM, beside the
## transient $U, one value per observation. TRE and TRE_Z have a persistent
## random effect but no persistent inefficiency, so they have no $H and asking
## for one is an error rather than a silent fall back to the transient score --
## which is exactly the component mix-up gap A9 shipped.
.te_persistent_slot <- function(object) object$H

## What this fit can actually answer, so the error never advises a `type` that
## has nothing to return either. The old message told a GTRE fit to use
## type = "jlms", where it stores no E[u | e] and fails again.
.eff_nothing_stored <- function(object, type, component) {
  mn <- if (is.null(object$model_name)) "unknown" else object$model_name
  if (identical(component, "persistent")) {
    return(paste0(
      "efficiency(component = \"persistent\"): model_name \"", mn, "\" stores no ",
      "persistent efficiency. Only the four-component panel models have one ",
      "(psfm() \"GTRE\", \"GTRE_FML\", \"GTRE_Z\"); \"TRE\" and \"TRE_Z\" have a ",
      "persistent random effect but no persistent INEFFICIENCY, and every ",
      "cross-sectional model has a single component. Use the default ",
      "component = \"transient\"."))
  }
  have <- c(
    if (!is.null(object$exp_u_hat) || !is.null(object$efficiency) ||
        !is.null(object$Efficiency) || !is.null(object$U)) "\"bc\"",
    if (!is.null(object$u_hat) || !is.null(object$jlms) ||
        (!is.null(object$u_posterior) && !is.null(object$u_posterior$mu_star))) "\"jlms\"",
    if (!is.null(object$u_posterior) && !is.null(object$u_posterior$mu_star)) "\"mode\"")
  paste0("efficiency(type = \"", type, "\"): model_name \"", mn, "\" stores ",
    if (identical(type, "bc")) "no E[exp(-u) | e]." else
      "neither E[u | e] nor the posterior needed to rebuild it.",
    if (length(have))
      paste0(" This fit can answer type = ", paste(have, collapse = " or "), ".")
    else " This fit stores no efficiency predictor at all.")
}

efficiency <- function(object,
                       type = c("bc", "jlms", "mode"),
                       logDepVar = TRUE,
                       newdata = NULL,
                       component = c("transient", "persistent")) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  type <- match.arg(type)
  .component_given <- !missing(component)
  component <- match.arg(component)
  if (length(logDepVar) != 1L || !is.logical(logDepVar) || is.na(logDepVar)) {
    stop("`logDepVar` must be TRUE or FALSE.", call. = FALSE)
  }

  ## The point estimate of u implied by each predictor. `bc` is a certainty
  ## equivalent rather than a moment of u, which is exactly why it differs from
  ## `jlms` by Jensen's inequality and not by noise.
  u <- switch(type,
    jlms = {
      ## Prefer the stored value; fall back to rebuilding it from the posterior
      ## with the SAME helper efficiency_ci() uses, so the two cannot drift.
      ## Most cross-sectional fits store exp_u_hat and the posterior but not
      ## E[u | e] itself.
      v <- object$u_hat
      if (is.null(v)) v <- object$jlms
      if (is.null(v)) {
        post <- object$u_posterior
        if (!is.null(post) && !is.null(post$mu_star)) {
          v <- .jlms_u(post$mu_star, post$sigma_star)
        }
      }
      if (is.null(v)) {
        ## NOT derived from the stored E[exp(-u) | e]. -log(E[exp(-u) | e]) is
        ## not E[u | e]; they differ by Jensen's inequality, and reporting one
        ## under the other's name is the defect gap A68 recorded for
        ## efficiency_ci(). A fit that stores only the Battese-Coelli score is
        ## told so, and told what it CAN answer.
        stop(.eff_nothing_stored(object, "jlms", component), call. = FALSE)
      }
      as.numeric(v)
    },
    bc = {
      v <- if (identical(component, "persistent")) .te_persistent_slot(object) else {
        ## Exhaustive rather than enumerated (gap A78). The chain used to end
        ## at $efficiency and so returned nothing for every psfm() GTRE/TRE
        ## entry, which store the components as $U and $H, nor for esfm(),
        ## which writes $Efficiency with a capital E. predict.sfareg() already
        ## carried the $U rung (R/sfareg_methods.R), so the two generics
        ## disagreed about whether a GTRE fit has an efficiency measure at all.
        w <- object$exp_u_hat
        if (is.null(w)) w <- object$efficiency
        if (is.null(w)) w <- object$Efficiency
        if (is.null(w)) w <- object$U
        w
      }
      if (is.null(v)) {
        stop(.eff_nothing_stored(object, "bc", component), call. = FALSE)
      }
      ## A four-component fit has TWO efficiencies, and "the" efficiency of
      ## such a fit is ambiguous. Said only when the caller did not choose --
      ## a user who passed component = explicitly already knows -- and only
      ## once a score has actually been found, so a call that is about to
      ## error does not also emit advice. Gap A9 was persistent inefficiency
      ## arriving inside a transient score, so silence here is the wrong
      ## default.
      if (!.component_given && identical(component, "transient") &&
          !is.null(.te_persistent_slot(object))) {
        message("efficiency(): model_name \"",
          if (is.null(object$model_name)) "unknown" else object$model_name,
          "\" reports two efficiencies. Returning the TRANSIENT one; ",
          "component = \"persistent\" gives exp(-h), one value per firm.")
      }
      -log(pmax(as.numeric(v), .SFA_CONSTANTS$MIN_POSITIVE))
    },
    mode = {
      post <- object$u_posterior
      if (is.null(post) || is.null(post$mu_star)) {
        stop("efficiency(type = \"mode\"): the mode needs the posterior of u ",
          "given the composed residual, which is a truncated normal only for ",
          "model_name \"NHN\", \"NHN_Z\", \"NE\" and \"NTN\". This fit is ",
          "model_name \"",
          if (is.null(object$model_name)) "unknown" else object$model_name,
          "\", whose posterior has a different shape.",
          call. = FALSE
        )
      }
      .u_mode(as.numeric(post$mu_star))
    }
  )

  if (isTRUE(logDepVar)) {
    ## y in logs: u is a proportional shortfall and TE = exp(-u). For `bc` this
    ## returns E[exp(-u) | e] unchanged, since the log and the exp cancel --
    ## deliberately, so the default path reports the stored number rather than
    ## a round-tripped copy of it.
    return(exp(-u))
  }

  ## A persistent score carries one value per FIRM and the fitted frontier one
  ## per OBSERVATION, so the level-scale ratio below has nothing to divide by.
  ## Recycling would misalign silently, which is how a firm-level quantity gets
  ## read as an observation-level one (see convergence/efficiency/, whose
  ## te_expand() refuses for the same reason).
  if (identical(component, "persistent")) {
    stop("efficiency(component = \"persistent\", logDepVar = FALSE): the ",
      "persistent score is one value per firm and the fitted frontier is one ",
      "per observation, so there is no level-scale ratio to form. Use ",
      "logDepVar = TRUE, which returns exp(-h) per firm.",
      call. = FALSE
    )
  }

  ## y in levels: the frontier is f(x) and TE = (f - u)/f = 1 - u/f. This is a
  ## RATIO of two estimated quantities, so it is less well behaved than the log
  ## version -- it is undefined where the fitted frontier is zero and negative
  ## where u exceeds it, both of which happen near zero output.
  ## .sfa_xb() rebuilds the frontier by re-evaluating the call's `data`
  ## argument, which fails whenever the fit was made inside a function whose
  ## environment has since gone -- a limitation shared by fitted(), residuals()
  ## and predict(). `newdata` is the escape hatch, and the error says so rather
  ## than reporting that the frontier does not exist.
  f <- tryCatch(as.numeric(.sfa_xb(object, newdata)$xb), error = function(e) NULL)
  if (is.null(f) || length(f) != length(u)) {
    stop("efficiency(logDepVar = FALSE): could not rebuild the fitted frontier ",
      "for this fit. This happens when the data used to fit it can no longer be ",
      "recovered from the call -- pass it explicitly as `newdata`.",
      call. = FALSE
    )
  }
  ## A cost frontier has y = f + u + v, so its efficiency is minimum over actual
  ## cost, f / (f + u) -- the level analogue of exp(-u) on the log scale. The
  ## production form 1 - u/f was applied to both until gap A41.
  cost <- .is_cost_fit(object)
  if (is.na(cost)) stop(.cost_fit_unknown("efficiency(logDepVar = FALSE)"), call. = FALSE)
  te <- if (cost) f / (f + u) else 1 - u / f
  if (any(!is.finite(te)) || any(te < 0, na.rm = TRUE)) {
    warning("efficiency(logDepVar = FALSE): ",
      sum(!is.finite(te) | te < 0, na.rm = TRUE),
      " of ", length(te), " scores are negative or non-finite, which happens ",
      "where the fitted frontier is at or below zero. Check that the dependent ",
      "variable really is on a level scale.",
      call. = FALSE
    )
  }
  te
}
