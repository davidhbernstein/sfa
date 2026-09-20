## S3 methods for class "sfareg"

coef.sfareg <- function(object, ...) {
  object$coefficients
}

vcov.sfareg <- function(object, type = c("hessian", "bhhh"), ...) {
  type <- match.arg(type)
  p <- length(object$coefficients)
  nm <- names(object$coefficients)

  ## BHHH: the inverse of the outer product of the per-observation scores.
  ## Worth having because it needs no Hessian at all -- it is defined whenever
  ## the scores are, which is exactly the case the default path fails on. When
  ## the Hessian is singular vcov() currently falls back to a DIAGONAL
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
        call. = FALSE
      )
    }
    dimnames(V) <- list(nm, nm)
    return(V)
  }

  if (!is.null(object$opt) && !is.null(object$opt$hessian)) {
    V <- tryCatch(solve(object$opt$hessian), error = function(e) NULL)
    if (!is.null(V)) {
      dimnames(V) <- list(nm, nm)
      return(V)
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

logLik.sfareg <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$value)) {
    warning("This fit has no stored optimizer output (e.g. psfm()'s GTRE_SEQ1/GTRE_SEQ2 are moment-based, not maximum likelihood), so logLik() is not defined for it.", call. = FALSE)
    ## Return a properly classed logLik carrying NA rather than a bare
    ## NA_real_.
    val <- NA_real_
    attr(val, "df") <- length(object$coefficients)
    attr(val, "nobs") <- nobs.sfareg(object)
    class(val) <- "logLik"
    return(val)
  }
  val <- -object$opt$value
  attr(val, "df") <- length(object$coefficients)
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
  ## The same identity problem .sfa_data() has (gap A47): the ordering below
  ## saves the function-local case, but a name that is absent from the
  ## formula's environment still falls through to whatever `dcall` happens to
  ## name in globalenv(). Validate what comes back before counting its rows,
  ## and report NA rather than a decoy's row count.
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

## Recovering the data a fit was built from (gap A47).
##
## The old route was eval(object$call$data, parent.frame(3)) -- a FIXED frame
## depth, correct only for the chain user -> fitted.sfareg -> .sfa_xb ->
## .sfa_data. Two things were wrong with it. It missed the fit's own frame
## whenever the model was fitted inside a function, so fitted(), residuals(),
## predict() and efficiency(logDepVar = FALSE) all failed on a perfectly good
## fit. Worse, when an UNRELATED object of the same name was visible at that
## depth it resolved to that object and returned it unchecked, so the methods
## answered silently from the wrong data. No fixed depth can be right for every
## call path, so the depth is no longer what is relied on: environments are
## searched in order of trust, and whatever comes back must identify itself as
## this fit's data before it is used.

## Environments a fit's `data` argument might resolve in, most trustworthy
## first. The formula's environment is where the model was specified, so it
## reaches the function-local frame a fixed parent.frame() depth cannot.
.sfa_data_envs <- function(object, frames = list()) {
  envs <- list()
  fenv <- if (inherits(object$formula, "formula")) environment(object$formula) else NULL
  if (is.environment(fenv)) envs <- c(envs, list(fenv))
  envs <- c(envs, Filter(is.environment, frames))
  c(envs, list(globalenv()))
}

## Does `dat` identify itself as the frame this fit was built from? Three
## tiers, strongest first; which are available depends on what the entry point
## stored. Returns the verdict and the name of the strongest check that could
## actually be applied, which the caller puts in its error message.
##
## Deliberately NOT used here: nobs.sfareg(). It re-evaluates the call itself
## when the fit carries no $nobs field, so validating against it would check a
## recovered frame against a count derived from the same recovery.
## The rows a fit actually USED, given a frame it was (or might have been)
## built from: complete cases across every pipe segment at once, response
## included, which is the rule data_proc2() applies at data.processing.R:219.
## Both the identity check and the design rebuild have to agree with it, or
## they disagree with the fit (gap A44).
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

  ## Tier 1, decisive wherever the fit kept its OLS residuals: rebuild them
  ## from the candidate and compare. sfm() stores them signed by `inefdec`
  ## (esfm.R), so either sign counts as a match. A tier that cannot be computed
  ## falls through to the weaker ones rather than rejecting -- an lm() that
  ## will not fit here is not evidence that the data is wrong.
  orr <- suppressWarnings(tryCatch(as.numeric(object$ols_residuals), error = function(e) NULL))
  if (length(orr)) {
    own <- tryCatch(
      {
        f1 <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
        ## On the rows the fit used, not on every row supplied: a row dropped
        ## only for a missing variance determinant is still complete on the
        ## frontier, so lm() would keep it and the lengths would disagree for
        ## the right data (gap A44).
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

  ## Nothing but structure was available. zsfm() and ttsfm() fits land here:
  ## they store no $nobs, no $data and no per-observation vector, so a
  ## same-shaped decoy cannot be told apart from the real thing. The trust
  ## ordering in .sfa_data_envs() is what protects those.
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

  ## A44: the fit dropped rows on complete.cases() across EVERY pipe segment at
  ## once, response included (data.processing.R:219). The rebuild below sees
  ## only the frontier, and delete.response() hides y from it, so a row missing
  ## a variance determinant -- or missing the response -- survives here although
  ## the fit never used it. The returned vector is then one row too long and
  ## silently misaligned against u_hat/exp_u_hat; residuals() went further and
  ## recycled the response against a shorter design. Apply the fit's own rule
  ## before building. Not applied to `newdata`, where the caller chooses the
  ## rows and there may be no response to be complete about.
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
