## influence_sfa(): which observations are driving this fit? Gap L7, from
## Stead, Wheat and Greene (2023), EJOR 309:188-201.
## See notes/code_history/influence.md.
##
## The package had no outlier or influence diagnostics at all, and the standard
## ones do not transfer: least-median-of-squares style outlier rules ignore the
## ASYMMETRY of the composed error, so a genuinely large u_i -- an inefficient
## firm, the thing being measured -- looks like an outlier and gets discarded.
## That defeats the purpose of the analysis.
##
## The right object is the influence function. For an M-estimator, and maximum
## likelihood is one, the influence function is a linear transformation of the
## score:
##
##   IF_i  proportional to  I(theta)^-1 psi(y_i, theta)
##
## and the estimator is B-robust exactly when that is BOUNDED. Both pieces are
## already computed here: `estfun.sfareg()` is the per-observation score matrix
## (added for I9) and `vcov()` is the inverse information. So this is assembly,
## as TIC() was.
##
## Stead, Wheat and Greene's substantive finding is that the canonical
## assumptions are NOT robust, while a Student's t noise term with fixed
## degrees of freedom satisfies their conditions for a wide range of
## inefficiency distributions. `sfm(model_name = "tHN")` is exactly that, so
## the package can already act on the diagnosis -- see the example in ?influence_sfa.

influence_sfa <- function(object, scale = TRUE) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  G <- estfun.sfareg(object)          # n x p per-observation scores
  ## G is on the ESTIMATION scale, so it has to be paired with a covariance on
  ## the SAME scale. vcov() reports on coef()'s scale, which for ttsfm(),
  ## copsfm() and ivsfm() is a different one -- and for ivsfm("IVLIML"), which
  ## estimates reduced-form parameters it never reports, a different LENGTH.
  ## Mixing the two silently rescales every influence value, and d_i below is
  ## only the metric-free quantity its comment claims when both agree. Which
  ## is why the covariance below is built from the Hessian directly for EVERY
  ## fit rather than only for the ones whose scales differ.
  ## The ESTIMATION-scale inverse observed information, UNMASKED, on both
  ## kinds of fit -- which is what bread() computes, so the two agree by
  ## construction.
  ##
  ## NOT stats::vcov(object). Its default path withholds a parameter's row and
  ## column wherever the fit reports no standard error (A77). That is right
  ## for REPORTED inference and wrong here: influence_sfa() is a diagnostic,
  ## and refusing to form it on the fits that most need diagnosing is the
  ## failure mode gap A81 records for plot(which = 1). Taking vcov() here made
  ## influence_sfa() error outright on psfm("GTRE", estimator = "sml") and on
  ## GTRE_FML, which it had always handled -- caught by CI on four platforms.
  ##
  ## .sfa_solve_equilibrated() rather than solve(): the two differ on exactly
  ## the fits it was written for -- GTRE_FML's Hessian carries a 2.5e17
  ## diagonal entry, so LAPACK refuses the direct solve while the scaled
  ## matrix is perfectly healthy -- so the old `else` branch failed there for
  ## a units artefact. Using one inverter on both branches also removes a
  ## difference between them that was never meant to be one.
  ## Still needed below, where IF is mapped onto the reported scale.
  jac <- object$par_scale
  idx <- object$par_index
  H <- if (!is.null(object$opt)) object$opt$hessian else NULL
  V <- if (is.null(H)) NULL else .sfa_solve_equilibrated(H)
  if (is.null(V) || any(!is.finite(V))) {
    stop("influence_sfa(): this fit's Hessian could not be inverted on the ",
      "estimation scale, so the influence function cannot be formed. Refit ",
      "with optHessian = TRUE, or try vcov(type = \"bhhh\"), which needs no ",
      "Hessian.",
      call. = FALSE
    )
  }
  n <- nrow(G)

  ## The empirical influence function. The n factor makes IF_i the
  ## approximate effect on the estimate of deleting observation i, which is
  ## what makes the numbers comparable to the coefficients themselves rather
  ## than merely rankable.
  IF <- G %*% V
  if (isTRUE(scale)) IF <- IF * n
  ## Map the influence onto the parameters the fit REPORTS -- the same delta
  ## method and marginalisation vcov() applies.
  if (is.matrix(jac)) {
    IF <- IF %*% t(jac)
  } else {
    if (!is.null(idx)) IF <- IF[, idx, drop = FALSE]
    if (!is.null(jac)) IF <- sweep(IF, 2L, jac, `*`)
  }
  colnames(IF) <- names(object$coefficients)

  ## A single number per observation. s_i' V s_i is the natural quadratic form
  ## here -- the same shape as Cook's distance, and metric-free, so it does not
  ## depend on how the parameters happen to be scaled.
  d <- rowSums((G %*% V) * G) * n

  ## Sensitivity: the sup-norm of the influence function. This is the quantity
  ## B-robustness is about. Two versions, and the distinction matters.
  ##
  ## The RAW sup-norm depends on how the parameters happen to be scaled, so it
  ## cannot be compared across models with different parameter vectors. Measured
  ## on n = 200: NHN 72.6, tHN 1322.5 -- which would say the robust
  ## specification is eighteen times worse, when all it reflects is tHN's extra
  ## nu on its own scale. Kept, because within one parameterisation (across
  ## tuning values, or clean against contaminated) it is the direct quantity.
  ##
  ## The SELF-STANDARDISED sup-norm, sup_i sqrt(d_i), measures the influence
  ## function in the information metric, so it is invariant to reparameterisation
  ## and IS comparable across models. On the same data, contaminating one
  ## response takes NHN from 13.5 to 32.3 and tHN from 13.1 to 17.2 -- which is
  ## Stead, Wheat and Greene's finding, and is invisible in the raw number.
  gamma_star <- max(sqrt(rowSums(IF^2)))
  gamma_std <- max(sqrt(d))
  per_par <- apply(abs(IF), 2, max)

  out <- list(
    influence = IF,
    d = d,
    sensitivity = gamma_star,
    sensitivity_std = gamma_std,
    max_abs_by_parameter = per_par,
    nobs = n,
    model_name = object$model_name,
    coefficients = object$coefficients,
    scaled = isTRUE(scale)
  )
  class(out) <- "sfa_influence"
  out
}

print.sfa_influence <- function(x, n_top = 5, ...) {
  cat("\nInfluence diagnostics for a stochastic frontier fit\n")
  cat("Stead, Wheat & Greene (2023), EJOR 309:188-201\n\n")
  cat("  model_name  : ", x$model_name, "\n", sep = "")
  cat("  observations: ", x$nobs, "\n", sep = "")
  cat(sprintf("  self-standardised sensitivity: %.4f\n", x$sensitivity_std))
  cat(sprintf("  raw sup |IF|                 : %.4f\n\n", x$sensitivity))

  ord <- order(-x$d)[seq_len(min(n_top, length(x$d)))]
  cat("Most influential observations:\n")
  d <- data.frame(obs = ord, d = round(x$d[ord], 4))
  print(d, row.names = FALSE)

  cat("\nLargest |influence| by parameter:\n")
  print(round(x$max_abs_by_parameter, 4))

  ## The interpretation that matters, and the one a user cannot get from the
  ## numbers alone: a large sensitivity is a property of the SPECIFICATION,
  ## not of the data. Deleting the offending observation does not fix a model
  ## whose influence function is unbounded -- the next most extreme point
  ## simply takes its place. The comparison to make is the self-standardised
  ## number BETWEEN specifications, not either number against a threshold:
  ## neither has a scale on which "large" means anything on its own.
  cat("\n  Sensitivity is a property of the specification, not of these data.\n")
  cat("  If it is large, deleting the top observations will not help: under\n")
  cat("  an unbounded influence function the next most extreme point takes\n")
  cat("  their place. Refit under a specification with bounded influence and\n")
  cat("  compare -- Stead, Wheat and Greene's remedy is a Student's t noise\n")
  cat("  term, model_name = \"tHN\" here; sfm(robust = ) is the other route.\n")
  cat("  Compare on the SELF-STANDARDISED number: the raw sup |IF| depends on\n")
  cat("  the parameterisation and is not comparable across models.\n\n")
  invisible(x)
}
