## S3 methods for class "sfareg"

coef.sfareg <- function(object, ...) {
  object$coefficients
}

## Map a covariance matrix from the ESTIMATION scale onto the scale coef()
## reports. Several entry points do not estimate what they report: ttsfm(),
## copsfm() and ivsfm() estimate LOG sigmas, ivsfm() reports
## rho = t / sqrt(1 + t't), and ivsfm("IVLIML") additionally estimates
## reduced-form parameters it never reports.
##
## Before this existed, vcov() named the estimation-scale matrix with the
## reported names and returned it: on a copsfm() fit sqrt(diag(vcov(f)))
## came back a factor of 1/sigma_u away from the sigma_u standard error the
## fit itself printed, and on ivsfm("IVLIML") -- 11 estimated, 6 reported --
## it failed outright in dimnames(). Both paths, Hessian and BHHH, go through
## here so the two can no longer disagree.
##
## Order matters. Subsetting is done on the INVERSE, because the marginal
## covariance of a sub-block is that block of the inverse, not the inverse of
## the block; the latter conditions on the nuisance parameters being known and
## understates every standard error.
.sfa_map_vcov <- function(V, idx, jac, p, nm, what) {
  if (is.matrix(jac)) {
    ## J is (reported x estimated), so it subsets and transforms at once --
    ## needed when a reported parameter is a non-diagonal function of several
    ## estimated ones, as ivsfm()'s rho is.
    if (ncol(jac) != nrow(V)) {
      stop(what, ": the stored parameter Jacobian is ", nrow(jac), "x",
        ncol(jac), " but the covariance has ", nrow(V),
        " parameters; refit under the current version.",
        call. = FALSE
      )
    }
    V <- jac %*% V %*% t(jac)
  } else {
    if (!is.null(idx)) {
      if (any(idx < 1L) || any(idx > nrow(V))) {
        stop(what, ": the stored parameter index does not match the ",
          "covariance; refit under the current version.",
          call. = FALSE
        )
      }
      V <- V[idx, idx, drop = FALSE]
    }
    if (!is.null(jac)) {
      if (length(jac) != nrow(V)) {
        stop(what, ": the stored parameter Jacobian does not match the ",
          "covariance; refit under the current version.",
          call. = FALSE
        )
      }
      ## Diagonal delta method: V * J J' elementwise is diag(J) V diag(J).
      V <- V * tcrossprod(jac)
    }
  }
  if (nrow(V) != p) {
    stop(what, ": the covariance maps onto ", nrow(V),
      " parameters but the fit reports ", p, ".",
      call. = FALSE
    )
  }
  dimnames(V) <- list(nm, nm)
  V
}

vcov.sfareg <- function(object, type = c("hessian", "bhhh", "sandwich", "clustered"),
                        cluster = NULL, ...) {
  type <- match.arg(type)
  p <- length(object$coefficients)
  nm <- names(object$coefficients)

  ## BHHH: the inverse of the outer product of the per-observation scores.
  ## Worth having because it needs no Hessian at all -- it is defined whenever
  ## the scores are, which is exactly the case the default path fails on. When
  ## the Hessian is singular vcov() otherwise falls back to a DIAGONAL
  ## approximation, discarding every covariance; BHHH keeps them.
  ##
  ## It is only as good as the information-matrix equality, so it disagrees
  ## with the Hessian under misspecification -- which is a feature when the two
  ## are compared deliberately and a trap when they are not. Not the default.
  if (identical(type, "bhhh")) {
    G <- tryCatch(estfun.sfareg(object), error = function(e) NULL)
    if (is.null(G)) {
      stop("vcov(type = \"bhhh\"): the score matrix is unavailable for this ",
        "fit. Refit with keep_objective = TRUE.",
        call. = FALSE
      )
    }
    V <- tryCatch(solve(crossprod(G)), error = function(e) NULL)
    if (is.null(V)) {
      stop("vcov(type = \"bhhh\"): the outer product of gradients is singular.",
        .sfa_dead_parameters(object, G),
        call. = FALSE
      )
    }
    return(.sfa_map_vcov(V, attr(G, "par_index"), attr(G, "par_scale"),
      p, nm, "vcov(type = \"bhhh\")"
    ))
  }

  ## Sandwich and clustered covariances, computed HERE rather than left to the
  ## sandwich package. Both are built on the ESTIMATION scale and only then
  ## mapped to the reported one, which is the step sandwich::sandwich() cannot
  ## take: it composes bread and meat itself and knows nothing about
  ## par_index/par_scale. Composing a reported-scale bread with an
  ## estimation-scale meat is what made those estimates silently wrong for
  ## every model whose two scales differ -- and impossible for ivsfm("IVLIML"),
  ## which estimates 11 parameters and reports 6, where the dimensions do not
  ## even conform.
  if (type %in% c("sandwich", "clustered")) {
    what <- sprintf("vcov(type = \"%s\")", type)
    G <- tryCatch(estfun.sfareg(object), error = function(e) NULL)
    if (is.null(G)) {
      stop(what, ": the score matrix is unavailable for this fit. Refit with ",
        "keep_objective = TRUE.",
        call. = FALSE
      )
    }
    ## Both failure paths below used to advise type = "bhhh" unconditionally,
    ## on the reasoning that it needs no Hessian. That is true and not
    ## sufficient: BHHH still has to invert crossprod(G), and the fits that
    ## reach these errors are frequently flat in a parameter, which kills the
    ## outer product too. Advising it blind sends the user to a second error.
    ## So ask first, and only recommend what is actually available.
    bhhh_ok <- !is.null(tryCatch(solve(crossprod(G)), error = function(e) NULL))
    bhhh_hint <- if (bhhh_ok) {
      " type = \"bhhh\" needs no Hessian and is defined here."
    } else {
      " type = \"bhhh\" needs no Hessian, but is ALSO undefined here: the outer product of the scores is singular too."
    }
    H <- if (!is.null(object$opt)) object$opt$hessian else NULL
    if (is.null(H)) {
      stop(what, ": this covariance needs the Hessian for its bread and this ",
        "fit carries none (was optHessian = FALSE?).", bhhh_hint,
        if (!bhhh_ok) .sfa_dead_parameters(object, G),
        call. = FALSE
      )
    }
    Hi <- .sfa_solve_equilibrated(H)
    if (is.null(Hi)) {
      stop(what, ": the Hessian is singular, so the bread is undefined.",
        bhhh_hint,
        .sfa_dead_parameters(object, G),
        call. = FALSE
      )
    }
    n <- nrow(G)
    if (identical(type, "sandwich")) {
      meat <- crossprod(G)
    } else {
      ## The unit a cluster indexes is the row of estfun(), which for psfm()
      ## and the other panel likelihoods is a FIRM, not a firm-year. Saying so
      ## is worth more than the length check on its own: an N*T column is the
      ## natural thing to reach for and it is the wrong length by construction.
      if (is.null(cluster)) {
        stop(what, ": `cluster` is required. It indexes the rows of the score ",
          "matrix, which is one row per ", if (!is.null(object$n_units)) {
            "FIRM"
          } else {
            "observation"
          }, " for this fit -- length ", n, ".",
          call. = FALSE
        )
      }
      cl <- if (is.data.frame(cluster)) interaction(cluster, drop = TRUE) else cluster
      cl <- as.factor(cl)
      if (length(cl) != n) {
        stop(what, ": `cluster` has length ", length(cl), " but the score ",
          "matrix has ", n, " rows. Scores are per ",
          if (!is.null(object$n_units)) "FIRM" else "observation",
          " for this fit, so the cluster vector must be that long and in the ",
          "same order.",
          call. = FALSE
        )
      }
      Gc <- rowsum(G, cl, reorder = FALSE)
      M <- nrow(Gc)
      if (M < 2L) {
        stop(what, ": `cluster` has a single group, so a clustered covariance ",
          "is not identified.",
          call. = FALSE
        )
      }
      ## Cluster adjustment M/(M-1) ONLY, which is what sandwich::vcovCL()
      ## applies by default -- it does not also apply the HC1 observation
      ## factor (n-1)/(n-k). Adding that factor made the two disagree by
      ## exactly (n-1)/(n-k), 1.4% on a 300-observation fit, which is how the
      ## difference was identified. Matching an established implementation is
      ## worth more here than picking the correction on first principles, and
      ## it lets the agreement be asserted in a test.
      meat <- crossprod(Gc) * (M / (M - 1))
    }
    V <- Hi %*% meat %*% Hi
    return(.sfa_map_vcov(V, attr(G, "par_index"), attr(G, "par_scale"),
      p, nm, what
    ))
  }

  if (!is.null(object$opt) && !is.null(object$opt$hessian)) {
    V <- .sfa_solve_equilibrated(object$opt$hessian)
    if (!is.null(V)) {
      V <- tryCatch(
        .sfa_map_vcov(V, object$par_index, object$par_scale, p, nm, "vcov()"),
        error = function(e) NULL
      )
      if (!is.null(V)) return(V)
    }
    warning("Hessian stored on this fit could not be inverted; falling back to a diagonal approximation built from the reported standard errors.", call. = FALSE)
  }

  if (!is.null(object$std.errors) && !all(is.na(object$std.errors))) {
    V <- diag(object$std.errors^2, nrow = p)
    dimnames(V) <- list(nm, nm)
    return(V)
  }

  warning("No Hessian or standard errors are available on this fit (was optHessian = FALSE?); returning a matrix of NAs.", call. = FALSE)
  matrix(NA_real_, p, p, dimnames = list(nm, nm))
}

## Degrees of freedom: the number of parameters ESTIMATED, which is not always
## the number reported. ivsfm("IVLIML") maximises over 11 and reports 6, so
## taking length(coefficients) charged AIC and BIC for 6 -- a free gift of
## 2 * 5 = 10 AIC points against any model that reports what it estimates.
## Falls back to the reported count for the estimators that carry no $opt.
## Names for the ESTIMATION-scale parameter vector, which is not always the
## reported one. ttsfm(), copsfm() and ivsfm() estimate log-sigmas and report
## sigmas; ivsfm("IVLIML") estimates 11 and reports 6. Anything indexed by
## opt$par -- the score matrix, the Hessian, a numerical gradient, a
## likelihood slice -- must be labelled with these, or entry j carries the
## name of a different parameter.
.sfa_est_names <- function(object) {
  par <- object$opt$par
  if (is.null(par)) {
    nm <- names(object$coefficients)
    return(if (is.null(nm)) {
      paste0("par", seq_along(object$coefficients))
    } else {
      nm
    })
  }
  p <- length(par)
  nm <- names(object$coefficients)
  if (length(nm) != p) {
    nm <- if (!is.null(names(par))) names(par) else paste0("par", seq_len(p))
  }
  nm
}

.sfa_npar <- function(object) {
  p <- if (!is.null(object$opt$par)) length(object$opt$par) else 0L
  if (p < 1L) length(object$coefficients) else as.integer(p)
}

logLik.sfareg <- function(object, ...) {
  ## The robust divergence estimators minimise an MLq / psi / density-power
  ## objective, NOT a negative log-likelihood, so -opt$value is not a log
  ## likelihood and AIC()/BIC() built on it are not comparable with anything.
  ## Measured on a 300-observation NHN fit: AIC came back -2471 against the
  ## MLE's +609, handing the robust fit a 3080-point advantage that is purely
  ## a change of objective. estfun() and TIC() already refuse these fits; this
  ## refuses them in the same terms rather than returning a plausible number.
  if (!is.null(object$robust) && !identical(object$robust, "mle")) {
    warning("logLik(): the robust divergence estimators (robust = ",
      dQuote(object$robust, FALSE), ") minimise a divergence, not a negative ",
      "log-likelihood, so logLik() -- and AIC()/BIC() with it -- is not ",
      "defined for this fit. Refit with robust = \"mle\" to compare models on ",
      "the likelihood.",
      call. = FALSE
    )
    val <- NA_real_
    attr(val, "df") <- .sfa_npar(object)
    attr(val, "nobs") <- nobs.sfareg(object)
    class(val) <- "logLik"
    return(val)
  }
  if (is.null(object$opt) || is.null(object$opt$value)) {
    warning("This fit has no stored optimizer output (e.g. psfm()'s GTRE_SEQ1/GTRE_SEQ2 are moment-based, not maximum likelihood), so logLik() is not defined for it.", call. = FALSE)
    ## Return a properly classed logLik carrying NA rather than a bare
    ## NA_real_.
    val <- NA_real_
    attr(val, "df") <- .sfa_npar(object)
    attr(val, "nobs") <- nobs.sfareg(object)
    class(val) <- "logLik"
    return(val)
  }
  ## opt$value is the OBJECTIVE the optimizer minimised, which is not always
  ## the negative log-likelihood. lcsfm("LCM_CN") with penalty_c > 0 maximises
  ## a penalised likelihood and stores the plain one separately; reporting the
  ## penalised objective here overstated logLik by the penalty and understated
  ## AIC by twice it -- a bias that always favoured the regularised fit, which
  ## is backwards, since a penalty is not evidence.
  val <- if (is.numeric(object$logLik_unpenalised) &&
    length(object$logLik_unpenalised) == 1L &&
    is.finite(object$logLik_unpenalised)) {
    object$logLik_unpenalised
  } else {
    -object$opt$value
  }
  attr(val, "df") <- .sfa_npar(object)
  attr(val, "nobs") <- nobs.sfareg(object)
  class(val) <- "logLik"
  val
}

nobs.sfareg <- function(object, ...) {
  ## Rows USED, not rows supplied. Re-evaluating the call (below) counts the
  ## latter, so a fit on data with any missing value reported too many
  ## observations and BIC() was computed against the wrong n.
  if (!is.null(object$nobs) && is.finite(object$nobs)) {
    return(as.integer(object$nobs))
  }
  if (!is.null(object$data)) {
    return(nrow(object$data))
  }
  ## Failing an explicit count, any vector the fit stores one element of per
  ## observation is exact where re-evaluating the call is not.
  for (nm in c("exp_u_hat", "u_hat", "residuals", "med_u_hat")) {
    v <- object[[nm]]
    if (is.numeric(v) && length(v) > 0L) {
      return(length(v))
    }
  }
  if (is.numeric(object$u_posterior$mu_star)) {
    return(length(object$u_posterior$mu_star))
  }
  ## sfm()/zsfm()/ttsfm() do not store the data on the fitted object, so fall
  ## back to re-evaluating the `data` argument of the recorded call.
  dcall <- object$call$data
  if (is.null(dcall)) {
    return(NA_integer_)
  }
  ## Validate what comes back before counting its rows: the ordering alone
  ## still lets a decoy in globalenv() through (gap A47).
  for (e in .sfa_data_envs(object, rev(sys.frames()))) {
    cand <- tryCatch(eval(dcall, envir = e), error = function(err) NULL)
    if (is.null(cand)) next
    if (!isTRUE(.sfa_data_check(object, cand)$ok)) next
    n <- tryCatch(nrow(as.data.frame(cand)), error = function(err) NULL)
    if (!is.null(n)) {
      return(as.integer(n))
    }
  }
  NA_integer_
}


## predict() / fitted() / residuals() for "sfareg".

## Environments a fit's `data` argument might resolve in, most trustworthy
## first: the formula's environment reaches a function-local fit, which a
## fixed parent.frame() depth cannot (gap A47).
.sfa_data_envs <- function(object, frames = list()) {
  envs <- list()
  fenv <- if (inherits(object$formula, "formula")) environment(object$formula) else NULL
  if (is.environment(fenv)) envs <- c(envs, list(fenv))
  envs <- c(envs, Filter(is.environment, frames))
  c(envs, list(globalenv()))
}

## The rows a fit actually used: complete cases across every pipe segment at
## once, response included, the rule data_proc2() applies (gap A44).
.sfa_complete_rows <- function(object, dat) {
  vars <- tryCatch(
    intersect(all.vars(stats::formula(Formula::Formula(object$formula))), names(dat)),
    error = function(e) character()
  )
  if (!length(vars)) {
    return(dat)
  }
  keep <- stats::complete.cases(dat[, vars, drop = FALSE])
  if (all(keep)) dat else dat[keep, , drop = FALSE]
}

## OLS residuals of the frontier on the rows a fit used, kept on zsfm() and
## ttsfm() fits as `anchor_resid` so .sfa_data_check() can rebuild them from
## recovered data and compare -- the decisive tier, where a row count alone
## passes any decoy of the same shape (gap A48). n doubles per fit.
.sfa_anchor_resid <- function(X, Y) {
  tryCatch(as.numeric(stats::lm.fit(as.matrix(X), as.numeric(Y))$residuals),
    error = function(e) NULL
  )
}

## Does `dat` identify itself as this fit's data? Three tiers, strongest
## first; which apply depends on what the entry point stored. nobs.sfareg() is
## deliberately not one of them -- it re-evaluates the call itself, so the
## check would be circular (gap A47).
.sfa_data_check <- function(object, dat) {
  no <- function(what) list(ok = FALSE, checked = what)
  dat <- tryCatch(as.data.frame(dat), error = function(e) NULL)
  if (is.null(dat) || !nrow(dat)) {
    return(no("none"))
  }

  ## Tier 3, always available: every variable the formula names is present.
  vars <- tryCatch(all.vars(stats::formula(Formula::Formula(object$formula))),
    error = function(e) character()
  )
  if (length(vars) && !all(vars %in% names(dat))) {
    return(no("variables"))
  }

  ## Tier 1, decisive where the fit kept its OLS residuals: rebuild and
  ## compare, either sign (sfm() stores them signed by `inefdec`). A tier that
  ## cannot be computed falls through rather than rejecting.
  ## zsfm() and ttsfm() keep the same residuals as `anchor_resid`: a private
  ## name, because skewness_test(), spec_test() and others read
  ## `ols_residuals` and would change behaviour on those fits (gap A48). Not
  ## `data_*`: `object$data` above partial-matches any such name.
  orr <- object$ols_residuals
  if (is.null(orr)) orr <- object$anchor_resid
  orr <- suppressWarnings(tryCatch(as.numeric(orr), error = function(e) NULL))
  if (length(orr)) {
    own <- tryCatch(
      {
        f1 <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
        ## On the rows the fit used: lm() would keep a row dropped only for a
        ## missing variance determinant, and the lengths would disagree.
        as.numeric(stats::resid(stats::lm(f1, data = .sfa_complete_rows(object, dat))))
      },
      error = function(e) NULL, warning = function(w) NULL
    )
    if (!is.null(own)) {
      if (length(own) != length(orr)) {
        return(no("ols_residuals"))
      }
      ok <- isTRUE(all.equal(own, orr, tolerance = 1e-6)) ||
        isTRUE(all.equal(-own, orr, tolerance = 1e-6))
      return(list(ok = ok, checked = "ols_residuals"))
    }
  }

  ## Tier 2: a count of rows USED, taken only from fields the fit stored
  ## outright.
  n_used <- NA_integer_
  if (!is.null(object$nobs) && length(object$nobs) == 1L && is.finite(object$nobs)) {
    n_used <- as.integer(object$nobs)
  }
  if (is.na(n_used)) {
    for (nm in c("exp_u_hat", "u_hat", "residuals", "med_u_hat")) {
      v <- object[[nm]]
      if (is.numeric(v) && length(v) > 0L) {
        n_used <- length(v)
        break
      }
    }
  }
  if (is.na(n_used) && is.numeric(object$u_posterior$mu_star)) {
    n_used <- length(object$u_posterior$mu_star)
  }
  if (!is.na(n_used)) {
    ## Rows SUPPLIED may exceed rows USED, because the fit dropped incomplete
    ## cases. A lower bound is the most this check can honestly assert.
    return(list(ok = nrow(dat) >= n_used, checked = "nobs"))
  }

  ## Structure only: a fit that stores none of the above. Every entry point's
  ## main path now stores at least `nobs`; a same-shaped decoy passes here.
  list(ok = TRUE, checked = "variables")
}

## Recover the data a fit was built from: explicit newdata wins, then the copy
## some models store on the object, then a validated search of the call.
.sfa_data <- function(object, newdata = NULL) {
  if (!is.null(newdata)) {
    return(as.data.frame(newdata))
  }
  if (!is.null(object$data)) {
    return(as.data.frame(object$data))
  }
  dcall <- object$call$data
  if (is.null(dcall)) {
    stop("Cannot recover the data used to fit this model. Pass `newdata` explicitly.",
      call. = FALSE
    )
  }
  frames <- rev(sys.frames())
  checked <- "none"
  for (e in .sfa_data_envs(object, frames)) {
    cand <- tryCatch(eval(dcall, envir = e), error = function(err) NULL)
    if (is.null(cand)) next
    chk <- .sfa_data_check(object, cand)
    checked <- chk$checked
    if (isTRUE(chk$ok)) {
      return(as.data.frame(cand))
    }
  }
  stop("Cannot recover the data used to fit this model: `", deparse(dcall),
    "` could not be found, or what was found is not the data this model was ",
    "fitted to (strongest available check: ", checked, "). Pass `newdata` explicitly.",
    call. = FALSE
  )
}

## Frontier design matrix and the matching beta, for any sfareg model.
.sfa_xb <- function(object, newdata = NULL) {
  dat <- .sfa_data(object, newdata)
  f1 <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)

  ## Apply the fit's own row rule before building, or the rebuild keeps rows
  ## the fit dropped and the result misaligns against u_hat (gap A44). Not for
  ## `newdata`: there the caller chooses the rows.
  if (is.null(newdata)) {
    dat <- .sfa_complete_rows(object, dat)
  }
  mm <- stats::model.matrix(stats::delete.response(stats::terms(f1, data = dat)), data = dat)
  cf <- object$coefficients
  nm <- intersect(colnames(mm), names(cf))
  if (!length(nm)) {
    stop("None of this fit's coefficients match the frontier design matrix; ",
      "cannot form a prediction.",
      call. = FALSE
    )
  }
  list(xb = as.numeric(mm[, nm, drop = FALSE] %*% cf[nm]), data = dat, formula = f1)
}

predict.sfareg <- function(object, newdata = NULL,
                           type = c("frontier", "response", "efficiency"), ...) {
  type <- match.arg(type)
  if (identical(type, "efficiency")) {
    if (!is.null(newdata)) {
      stop("type = \"efficiency\" is only available for the estimation sample: ",
        "predicting efficiency requires the composed residual, which needs ",
        "the response.",
        call. = FALSE
      )
    }
    te <- object$exp_u_hat
    if (is.null(te)) te <- object$U
    if (is.null(te)) {
      stop("This model (", object$model_name, ") does not return an efficiency ",
        "prediction; see ?sfm for which models do.",
        call. = FALSE
      )
    }
    return(as.numeric(te))
  }
  z <- .sfa_xb(object, newdata)
  if (identical(type, "frontier")) {
    return(z$xb)
  }

  ## "response": the frontier shifted by predicted inefficiency, i.e.
  u <- object$u_hat
  if (is.null(u) && !is.null(object$exp_u_hat)) u <- -log(pmax(object$exp_u_hat, .Machine$double.xmin))
  if (is.null(u)) {
    stop("type = \"response\" needs an inefficiency prediction, which model ",
      object$model_name, " does not return.",
      call. = FALSE
    )
  }
  if (!is.null(newdata)) {
    stop("type = \"response\" is only available for the estimation sample.", call. = FALSE)
  }
  cost <- .is_cost_fit(object)
  if (is.na(cost)) stop(.cost_fit_unknown("predict(type = \"response\")"), call. = FALSE)
  z$xb - (if (cost) -1 else 1) * as.numeric(u)
}

## TRUE for a cost frontier (inefdec = FALSE), FALSE for production (the
## default), NA if the call's `inefdec` is a name that can no longer be
## evaluated. It is looked up where the formula was created, which is where a
## variable passed as `inefdec` normally lives.
.is_cost_fit <- function(object) {
  v <- object$call$inefdec
  if (is.null(v)) {
    return(FALSE)
  }
  if (!is.logical(v)) {
    env <- tryCatch(environment(stats::formula(object$formula)), error = function(e) NULL)
    v <- tryCatch(eval(v, if (is.null(env)) globalenv() else env), error = function(e) NA)
  }
  if (length(v) != 1L || !is.logical(v) || is.na(v)) {
    return(NA)
  }
  !v
}

.cost_fit_unknown <- function(what) {
  paste0(what, ": cannot tell whether this fit is a production or a cost ",
    "frontier, because its call gives `inefdec` as a name that no longer ",
    "evaluates. Refit with inefdec = TRUE or FALSE written out.")
}

fitted.sfareg <- function(object, ...) .sfa_xb(object)$xb

residuals.sfareg <- function(object, ...) {
  z <- .sfa_xb(object)
  yv <- all.vars(z$formula)[1]
  if (!(yv %in% names(z$data))) {
    stop("Response '", yv, "' not found in the fit's data; cannot form residuals.",
      call. = FALSE
    )
  }
  as.numeric(z$data[[yv]]) - z$xb
}
