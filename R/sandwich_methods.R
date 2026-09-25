## bread() and estfun() for "sfareg", which is all the `sandwich` package needs
## to work against these fits: vcovHC() for heteroskedasticity-consistent
## standard errors, vcovCL() for clustered ones, and coeftest() to print them.
## See notes/code_history/sandwich_methods.md.
##
## vcov.sfareg() alone is valid only under correct specification and independence.

## The "bread" of the sandwich: n times the inverse of the negative Hessian of
## the summed log-likelihood, i.e. n * vcov.
bread.sfareg <- function(x, ...) {
  V <- stats::vcov(x)
  ## psfm() scores are per firm, so the bread must scale by the same unit
  ## count estfun() returns rows for; nobs() (N*T) would inflate it by T^2.
  n <- if (!is.null(x$n_units)) x$n_units else stats::nobs(x)
  if (!is.finite(n)) {
    stop("bread(): the number of observations is unavailable for this fit.",
      call. = FALSE
    )
  }
  V * n
}

## Per-observation scores d loglik_i / d theta (n x p), by central differences of the
## likelihood kept by keep_objective = TRUE, step scaled per parameter.
estfun.sfareg <- function(x, ...) {
  if (is.null(x$objective)) {
    stop("estfun(): this fit does not retain its likelihood, so the score ",
      "matrix cannot be built. Refit with keep_objective = TRUE: that works for ",
      "sfm() fits, and for psfm() only with \"GTRE\" (estimator = \"sml\"), ",
      "\"TRE\", \"GTRE_Z\" and \"TRE_Z\". Other psfm() models and the other ",
      "entry points have no per-observation likelihood.",
      call. = FALSE
    )
  }
  if (!is.null(x$robust) && !identical(x$robust, "mle")) {
    stop("estfun(): the robust divergence estimators (robust = ",
      dQuote(x$robust), ") do not maximise a log-likelihood, so a score ",
      "matrix is not defined for them. sandwich-based standard errors do not ",
      "apply; sfm() already reports sandwich-form errors for those fits.",
      call. = FALSE
    )
  }
  par <- x$opt$par
  if (is.null(par)) {
    stop("estfun(): this fit has no parameter vector (was it produced by a ",
      "non-likelihood estimator?).",
      call. = FALSE
    )
  }
  ## Column names come from the ESTIMATION-scale vector, which is not always
  ## the reported one. ttsfm(), copsfm() and ivsfm() estimate log-sigmas and
  ## report sigmas, and ivsfm("IVLIML") estimates reduced-form nuisance
  ## parameters it never reports -- so names(coefficients) can be the wrong
  ## length or the wrong scale. Falling back to it blindly is how a score
  ## matrix on the log scale gets labelled as if it were on the level scale.
  p <- length(par)
  nm <- names(x$coefficients)
  if (length(nm) != p) {
    nm <- if (!is.null(names(par))) names(par) else paste0("par", seq_len(p))
  }

  ll_i <- function(theta) {
    v <- tryCatch(x$objective(theta, per_obs = TRUE), error = function(e) NULL)
    if (is.null(v)) {
      ## psfm()'s closures were never given the `per_obs` branch that sfm()'s
      ## gained in 1.2.0, so keep_objective = TRUE stores an objective the
      ## score matrix cannot use. Say which case this is rather than blaming a
      ## stale fit for both.
      stop("estfun(): the stored likelihood does not return per-observation ",
        "contributions, so the score matrix cannot be built.\n",
        if (!is.null(x$model_name)) {
          paste0("  This fit is model_name = ", dQuote(x$model_name), ". ")
        } else {
          "  "
        },
        "Per-observation contributions are available for sfm() fits and for ",
        "psfm()'s GTRE, TRE, GTRE_Z and TRE_Z (per firm) from 1.2.1. A fit ",
        "retained by an older version must be refit.",
        call. = FALSE
      )
    }
    as.numeric(v)
  }

  base <- ll_i(par)
  n <- length(base)
  eps <- .Machine$double.eps^(1 / 3)
  sc <- matrix(NA_real_, n, p, dimnames = list(NULL, nm))
  for (j in seq_len(p)) {
    h <- eps * max(abs(par[j]), 1)
    up <- dn <- par
    up[j] <- par[j] + h
    dn[j] <- par[j] - h
    sc[, j] <- (ll_i(up) - ll_i(dn)) / (2 * h)
  }
  ## A non-finite score would silently poison the whole meat matrix.
  sc[!is.finite(sc)] <- 0
  ## How the reported parameters sit inside this estimation-scale vector.
  ## vcov(type = "bhhh") needs both to get from here to a covariance on the
  ## scale coef() reports: invert on THIS scale, then subset, then delta.
  if (!is.null(x$par_index)) attr(sc, "par_index") <- x$par_index
  if (!is.null(x$par_scale)) attr(sc, "par_scale") <- x$par_scale
  sc
}
