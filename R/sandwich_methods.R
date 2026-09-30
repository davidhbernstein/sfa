## bread() and estfun() for "sfareg", which is all the `sandwich` package needs
## to work against these fits: vcovHC() for heteroskedasticity-consistent
## standard errors, vcovCL() for clustered ones, and coeftest() to print them.
## See notes/code_history/sandwich_methods.md.
##
## vcov.sfareg() alone is valid only under correct specification and independence.

## The "bread" of the sandwich: n times the inverse of the negative Hessian of
## the summed log-likelihood, i.e. n * vcov.
bread.sfareg <- function(x, ...) {
  ## ON THE ESTIMATION SCALE, to match estfun(). This used to return
  ## vcov(x) * n, which is the REPORTED scale, and sandwich::sandwich()
  ## composes the two into bread %*% meat %*% bread without knowing they
  ## disagree. Where the two scales have the same dimension that produced
  ## finite, plausible, WRONG numbers -- on ttsfm("TTNE") the sigma standard
  ## errors came out 29% low, 18% low and 26% high, each in a different
  ## direction, while the betas (whose Jacobian is 1) were right. Where the
  ## dimensions differ it produced "non-conformable arguments", which is how
  ## ivsfm("IVLIML") surfaced it.
  ##
  ## So sandwich::sandwich() and sandwich::vcovCL() now return an
  ## ESTIMATION-scale covariance: internally consistent, and for IVLIML
  ## defined at all. For a covariance on the scale coef() reports, use
  ## vcov(type = "sandwich") or vcov(type = "clustered"), which map through
  ## par_index and par_scale afterwards. For every model whose two scales
  ## agree -- sfm(), zsfm(), lcsfm(), psfm()'s TRE/GTRE family -- the two are
  ## identical and nothing changes.
  n <- if (!is.null(x$n_units)) x$n_units else stats::nobs(x)
  if (!is.finite(n)) {
    stop("bread(): the number of observations is unavailable for this fit.",
      call. = FALSE
    )
  }
  H <- if (!is.null(x$opt)) x$opt$hessian else NULL
  if (!is.null(H)) {
    Hi <- .sfa_solve_equilibrated(H)
    if (!is.null(Hi)) {
      nmE <- .sfa_est_names(x)
      dimnames(Hi) <- list(nmE, nmE)
      return(Hi * n)
    }
  }
  ## No usable Hessian. Falling back to the reported-scale covariance would
  ## reintroduce the mismatch silently, so only do it where the two scales are
  ## known to agree.
  if (!is.null(x$par_index) || !is.null(x$par_scale)) {
    stop("bread(): this fit has no invertible Hessian, and its estimation and ",
      "reported parameter scales differ, so no estimation-scale bread can be ",
      "formed. Use vcov(type = \"bhhh\"), which needs no Hessian.",
      call. = FALSE
    )
  }
  stats::vcov(x) * n
}

## Invert a Hessian after equilibrating it to unit diagonal.
##
## solve(H) fails on these fits for a reason that is arithmetic, not
## statistical. A variance component pinned near its floor makes the likelihood
## enormously sharp in that one direction: on psfm("GTRE_FML") at N = 40, T = 5,
## diag(H) runs 8e1 ... 1e3 for six parameters and 2.5e17 for the seventh, so
## the reciprocal condition is 6e-17 and LAPACK refuses. The SCALED matrix is
## perfectly healthy -- eigenvalues 0.13 to 2.4 -- so there is no flat
## direction and nothing is unidentified. It is a units artefact.
##
## H^-1 = D^-1/2 (D^-1/2 H D^-1/2)^-1 D^-1/2 is exact, so this changes no
## answer where the direct solve already worked: measured agreement 2e-14 on
## seeds where both succeed. Where the direct solve failed it returns finite
## standard errors for every parameter instead of a row of NaN, with the
## boundary parameter's own error correctly near zero.
##
## Falls back to the plain solve when the diagonal is not usable (a
## non-positive curvature is a real problem, not a scaling one).
.sfa_solve_equilibrated <- function(H) {
  if (!is.matrix(H) || nrow(H) != ncol(H) || !all(is.finite(H))) return(NULL)
  dg <- diag(H)
  if (all(is.finite(dg)) && all(dg > 0)) {
    sc <- 1 / sqrt(dg)
    Vs <- tryCatch(solve(H * outer(sc, sc)), error = function(e) NULL)
    if (!is.null(Vs)) {
      V <- Vs * outer(sc, sc)
      if (all(is.finite(V))) return(V)
    }
  }
  tryCatch(solve(H), error = function(e) NULL)
}

## Name the parameters a singular covariance is singular BECAUSE of.
##
## "the outer product of gradients is singular" is true and useless. When a
## parameter is unidentified at the fitted point the likelihood is flat in it,
## its score column is numerically dead, and NO covariance estimator can supply
## a standard error for it -- so the useful thing to report is WHICH parameter,
## not that a matrix inversion failed. The evidence is cheap: a score column at
## 1e-9 beside neighbours at 1e1, a column at 1e15 that will cancel
## catastrophically in crossprod(), or a coefficient pinned at machine epsilon.
##
## Observed in three places before this was written: psfm("GTRE_FML") with
## sigr at 2.22e-16 and two score columns exactly zero beside two at 1e15;
## psfm("K1990") with the decay parameters saturated at b = -17.3, c = -4.0 and
## scores ~2.5e-9 against ~10-35; and sfm("NNAK") returning one constant
## efficiency for every firm. See horserace/FUNCTIONALITY_GAPS.md A55.
## H, when supplied, is read for its near-null eigenvectors. That is the
## name-free way to answer "which parameter died": a flat direction IS a small
## eigenvalue, and the loadings of its eigenvector say which coefficients move
## along it. Reading it off the coefficients instead would mean guessing which
## ones are scales from their names, which is fragile and wrong for the models
## that report a lambda/sigma reparameterization.
.sfa_dead_parameters <- function(x, G = NULL, H = NULL) {
  nm <- tryCatch(.sfa_est_names(x), error = function(e) NULL)
  dead <- character(0)
  huge <- character(0)
  pinned <- character(0)
  if (!is.null(G) && is.matrix(G) && ncol(G) > 0L) {
    cmax <- apply(abs(G), 2L, function(z) max(z[is.finite(z)], na.rm = TRUE))
    cmax[!is.finite(cmax)] <- 0
    ref <- stats::median(cmax[cmax > 0])
    if (is.finite(ref) && ref > 0) {
      lab <- if (!is.null(nm) && length(nm) == length(cmax)) nm else
        paste0("[", seq_along(cmax), "]")
      dead <- lab[cmax < 1e-6 * ref]
      huge <- lab[cmax > 1e6 * ref]
    }
  }
  par <- if (!is.null(x$opt)) x$opt$par else NULL
  if (!is.null(par) && length(par)) {
    lab <- if (!is.null(nm) && length(nm) == length(par)) nm else
      paste0("[", seq_along(par), "]")
    pinned <- lab[is.finite(par) & abs(par) <= 1e-10]
  }
  msg <- character(0)
  if (length(dead)) {
    msg <- c(msg, sprintf("the likelihood is flat in %s (score column%s numerically zero)",
      paste(dQuote(dead, FALSE), collapse = ", "), if (length(dead) > 1L) "s" else ""))
  }
  if (length(pinned)) {
    msg <- c(msg, sprintf("%s sits at zero to machine precision",
      paste(dQuote(pinned, FALSE), collapse = ", ")))
  }
  if (length(huge)) {
    msg <- c(msg, sprintf("%s carries a score of order %.0e, which cancels catastrophically in the crossproduct",
      paste(dQuote(huge, FALSE), collapse = ", "),
      max(apply(abs(G[, match(huge, if (!is.null(nm) && length(nm) == ncol(G)) nm else paste0("[", seq_len(ncol(G)), "]")), drop = FALSE]), 2L, max))))
  }
  if (!is.null(H) && is.matrix(H) && nrow(H) == ncol(H) && all(is.finite(H))) {
    ## NORMALISE BY THE DIAGONAL FIRST. The raw Hessian here is badly scaled:
    ## a variance component pinned near zero makes the likelihood enormously
    ## SHARP in that direction, so its eigenvalue dwarfs the rest and a
    ## threshold against max(eigenvalue) then calls every OTHER direction flat.
    ## An earlier version of this did exactly that and named all six healthy
    ## parameters while omitting the one that had collapsed -- worse than
    ## silence, because it reads like an answer. Scaling to unit diagonal makes
    ## the eigenvalues comparable across parameters, so a small one means a
    ## genuinely flat direction rather than a units artefact.
    dg <- diag(H)
    lab <- if (!is.null(nm) && length(nm) == nrow(H)) nm else
      paste0("[", seq_len(nrow(H)), "]")
    if (all(is.finite(dg)) && all(dg > 0)) {
      sc <- 1 / sqrt(dg)
      Hs <- H * outer(sc, sc)
      ev <- tryCatch(eigen(Hs, symmetric = TRUE), error = function(e) NULL)
      if (!is.null(ev)) {
        av <- abs(ev$values)
        flat <- which(av < 1e-8)
        if (length(flat)) {
          who <- unique(unlist(lapply(flat, function(k) {
            w <- abs(ev$vectors[, k]); lab[w >= 0.5 * max(w)]
          })))
          msg <- c(msg, sprintf(
            "the likelihood is flat in %d direction%s of the scaled Hessian, loading on %s",
            length(flat), if (length(flat) > 1L) "s" else "",
            paste(dQuote(who, FALSE), collapse = ", ")))
        }
      }
    }
    ## A non-positive curvature on the diagonal is its own diagnosis: the fit
    ## is not at a maximum in that parameter.
    bad <- lab[is.finite(dg) & dg <= 0]
    if (length(bad)) {
      msg <- c(msg, sprintf("the Hessian has non-positive curvature in %s, so the fit is not at a maximum in %s",
        paste(dQuote(bad, FALSE), collapse = ", "),
        if (length(bad) > 1L) "those directions" else "that direction"))
    }
  }
  if (!length(msg)) return(NULL)
  paste0(" Likely cause: ", paste(msg, collapse = "; "),
    ". A parameter the likelihood is flat in has no standard error under any ",
    "estimator; refit without it, or on a design that identifies it.")
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
