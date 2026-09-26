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
  if (!is.null(x$robust) && !identical(x$robust, "mle")) {
    stop("estfun(): the robust divergence estimators (robust = ",
      dQuote(x$robust), ") do not maximise a log-likelihood, so a score ",
      "matrix is not defined for them. sandwich-based standard errors do not ",
      "apply; sfm() already reports sandwich-form errors for those fits.",
      call. = FALSE
    )
  }
  if (is.null(x$objective)) {
    stop("estfun(): this fit does not retain its likelihood, so the score ",
      "matrix cannot be built. Refit with keep_objective = TRUE. That works ",
      "for every maximum-likelihood fit: sfm() (except its robust ",
      "divergence estimators), zsfm(), ",
      "lcsfm(), ttsfm() (\"TTNE\", \"TTHN\"), copsfm(), selsfm(), ivsfm() ",
      "(\"IVLIML\", \"IVCF\"), and psfm() with ",
      paste(dQuote(.PSFM_SCORE_MODELS, FALSE), collapse = ", "),
      ". The remaining psfm() models, ttsfm(\"TTNLS\") and ivsfm(\"C2SLS\") ",
      "are not fit by maximum likelihood and have no score matrix.",
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
  nm <- .sfa_est_names(x)

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
        "Every maximum-likelihood entry point supplies them as of 1.2.1, and ",
        "psfm() supplies them per FIRM. A fit retained by an older version ",
        "must be refit under the current one.",
        call. = FALSE
      )
    }
    as.numeric(v)
  }

  base <- ll_i(par)
  n <- length(base)
  eps <- .Machine$double.eps^(1 / 3)
  sc <- matrix(NA_real_, n, p, dimnames = list(NULL, nm))
  ## A bail penalty is a SENTINEL, not a likelihood value. Every closure
  ## returns one when a step leaves the admissible region, and the smallest in
  ## use is MAX_VALUE^0.1 ~ 6e30 -- finite, so it passes is.finite() and
  ## differences to a score of order 1e150 that looks like a number. Catch it
  ## by magnitude.
  ##
  ## The threshold is that penalty divided by n, because the per-observation
  ## form SPREADS it: copsfm() and ivsfm() return rep(-PEN / n, n) while
  ## psfm() returns rep(-PEN, n) undivided. Testing against PEN itself misses
  ## the spread form entirely -- which it did, silently, until a synthetic
  ## likelihood in test-sandwich.R was written to force the branch. Even
  ## divided this sits ~25 orders above any real contribution.
  pen_floor <- .SFA_CONSTANTS$MAX_VALUE^0.1 / max(n, 1L)
  n_bailed <- 0L
  for (j in seq_len(p)) {
    h <- eps * max(abs(par[j]), 1)
    up <- dn <- par
    up[j] <- par[j] + h
    dn[j] <- par[j] - h
    vu <- ll_i(up)
    vd <- ll_i(dn)
    bu <- !is.finite(vu) | abs(vu) >= pen_floor
    bd <- !is.finite(vd) | abs(vd) >= pen_floor
    d <- (vu - vd) / (2 * h)
    ## One side outside the region: difference against the centre instead of
    ## against the penalty, which is the same derivative to first order.
    fwd <- bd & !bu
    bwd <- bu & !bd
    if (any(fwd)) d[fwd] <- (vu[fwd] - base[fwd]) / h
    if (any(bwd)) d[bwd] <- (base[bwd] - vd[bwd]) / h
    ## Both sides outside: the parameter is at an edge and no derivative
    ## exists. Zero is the same choice the non-finite guard below makes.
    d[bu & bd] <- 0
    n_bailed <- n_bailed + sum(bu | bd)
    sc[, j] <- d
  }
  ## A non-finite score would silently poison the whole meat matrix.
  sc[!is.finite(sc)] <- 0
  ## Recorded so a caller can tell "the scores vanish at the optimum" from
  ## "the optimum is on a boundary, where they need not". The first-order
  ## conditions do not hold at a bound, so a test of them has to know.
  attr(sc, "n_bailed") <- n_bailed
  if (n_bailed > 0L) {
    warning("estfun(): ", n_bailed, " finite-difference evaluation",
      if (n_bailed == 1L) "" else "s",
      " left the admissible region, so this fit sits at or near a parameter ",
      "bound. Those scores use a one-sided difference, or are zero where ",
      "both sides are inadmissible; OPG/BHHH standard errors built from them ",
      "understate the uncertainty in the affected directions.",
      call. = FALSE
    )
  }
  ## How the reported parameters sit inside this estimation-scale vector.
  ## vcov(type = "bhhh") needs both to get from here to a covariance on the
  ## scale coef() reports: invert on THIS scale, then subset, then delta.
  if (!is.null(x$par_index)) attr(sc, "par_index") <- x$par_index
  if (!is.null(x$par_scale)) attr(sc, "par_scale") <- x$par_scale
  sc
}
