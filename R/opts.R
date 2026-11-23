## Convergence code for "the final optim() polish stage did not run or did not
## produce a usable result". optim() itself returns 0, 1, 10, 51 or 52, so 99
## cannot collide with one of its own codes.
.SFA_OPTIM_SKIPPED <- 99L

## Does this point look like a real optimum, or like a boundary collapse? Used
## only by the A42 guard in opt.optim(); it must reject the NR spike that
## optim() is right to escape and accept an ordinary interior point. See
## notes/code_history/opts.md.
.sfa_point_credible <- function(fn, par, val, lower = NA, upper = NA) {
  if (!is.finite(val)) return(FALSE)
  ## Exactly 0 catches -0. A summed log-likelihood of exactly zero means every
  ## density evaluated to exactly 1: a collapse, not a fit (gap A49).
  if (val == 0) return(FALSE)
  lo <- rep_len(lower, length(par))
  up <- rep_len(upper, length(par))
  ## AT a bound, or BEYOND it. The two used to be different: these tests read
  ## `abs(par - lo) <= tol`, which rejects a point resting on the floor and
  ## accepts one sitting below it -- so the gate refused the boundary and
  ## admitted the region past it, which is strictly worse. psfm("TRE_Z") is
  ## how that surfaced: sigma_r = -0.2 against a floor of MIN_POSITIVE is
  ## 0.2 away from it, far outside the tolerance, so a negative standard
  ## deviation passed as a credible point. `<=` in place of `abs(...) <=`
  ## covers both cases and is the whole fix.
  tol_lo <- 1e-8 * pmax(1, abs(lo))
  tol_up <- 1e-8 * pmax(1, abs(up))
  if (any(is.finite(lo) & par <= lo + tol_lo)) return(FALSE)
  if (any(is.finite(up) & par >= up - tol_up)) return(FALSE)
  pc <- tryCatch(suppressWarnings(fn(par, per_obs = TRUE)), error = function(e) NULL)
  if (!is.null(pc) && length(pc) > 1L && all(is.finite(pc))) {
    if (all(pc == 0)) return(FALSE)
    if (min(pc) == max(pc)) return(FALSE)
  }
  TRUE
}

## Put a stage-1/2 optimizer result into optim()'s shape before it becomes a
## fit's `$opt` (gap A52). Every consumer of `$opt` -- logLik(), print(),
## summary(), sfa_diagnostics(), lcsfm()'s penalty arithmetic -- reads
## optim()'s field names, and two of the four stages do not use them:
##
##   optim    par value counts convergence message hessian   reference
##   psoptim  par value counts convergence message           conforms
##   nlminb   par objective convergence iterations evaluations message
##   bobyqa   par fval feval ierr msg
##
## Stored raw, a bobyqa object made logLik() NA and a lcsfm() fit a hard error
## (`-opt$value` on NULL). Conforming objects pass through unchanged, so the
## optHessian = TRUE path -- the default -- is untouched. There is no hessian
## here: these are the optHessian = FALSE paths, which report none.
.as_optim <- function(z) {
  if (is.null(z) || !is.null(z$value)) {
    return(z)
  }
  if (!is.null(z$objective)) {
    ev <- z$evaluations
    return(list(
      par = z$par,
      value = z$objective,
      counts = c(`function` = unname(ev[1]), gradient = unname(ev[2])),
      convergence = as.integer(z$convergence),
      message = z$message
    ))
  }
  if (!is.null(z$fval)) {
    return(list(
      par = z$par,
      value = z$fval,
      counts = c(`function` = as.integer(z$feval), gradient = NA_integer_),
      convergence = as.integer(z$ierr),
      message = z$msg
    ))
  }
  z
}

## opt.nlminb() -- primary optimization stage.
opt.nlminb <- function(fn, start_v, lower.nlminb, upper.nlminb = Inf,
                       gr = NULL, maxit.nlminb = 500, nlminb.TF = TRUE,
                       verbose = FALSE) {
  start_feval <- fn(start_v)
  nlm1 <- NULL
  if (isTRUE(nlminb.TF)) {
    nlm1 <- tryCatch(
      stats::nlminb(
        start = start_v, objective = fn, gradient = gr,
        lower = lower.nlminb, upper = upper.nlminb,
        control = list(
          iter.max = maxit.nlminb,
          eval.max = 2L * maxit.nlminb,
          trace = if (verbose) 1L else 0L
        )
      ),
      error = function(e) NULL
    )
    ## Only accept the new point if it actually improved the objective.
    if (!is.null(nlm1) && is.finite(nlm1$objective) && isTRUE(start_feval > nlm1$objective)) {
      start_v <- nlm1$par
      start_feval <- nlm1$objective
    }
  }
  results <- list(start_v, start_feval, nlm1)
  names(results) <- c("start_v", "start_feval", "nlm1")
  return(results)
}

opt.bobyqa <- function(fn, start_v, lower.bobyqa, upper.bobyqa = Inf, maxit.bobyqa, bob.TF, rhobeg = NA, rhoend = NA, verbose = verbose) {
  if (isTRUE(bob.TF == TRUE)) {
    ## bobyqa() stops on a start outside its bounds, and a boundary
    ## random-effects start can put a scale below its 1e-7 bound on some
    ## platforms (gap A35). Only coordinates strictly outside are moved.
    lo_v <- rep_len(lower.bobyqa, length(start_v))
    up_v <- rep_len(upper.bobyqa, length(start_v))
    below <- is.finite(lo_v) & start_v < lo_v
    above <- is.finite(up_v) & start_v > up_v
    start_v[below] <- lo_v[below] + 1e-6 * pmax(1, abs(lo_v[below]))
    start_v[above] <- up_v[above] - 1e-6 * pmax(1, abs(up_v[above]))
  }
  start_feval <- fn(start_v)
  bob1 <- NULL
  if (isTRUE(bob.TF == TRUE)) {
    bob1 <- bobyqa(
      par = start_v,
      fn = fn,
      lower = lower.bobyqa,
      upper = upper.bobyqa,
      control = list(
        iprint = if (verbose) 2 else 0,
        maxfun = maxit.bobyqa,
        rhobeg = rhobeg,
        rhoend = rhoend
      )
    )

    if (isTRUE(start_feval > bob1$fval)) {
      start_v <- bob1$par
      start_feval <- fn(start_v)
    }
  }

  results <- list(start_v, start_feval, bob1)
  names(results) <- c("start_v", "start_feval", "bob1")
  return(results)
}

## Two points closer than this (per coordinate, relative to max(1, |x|)) share
## a Hessian: 1e-8 is 1e-5 of optim()'s default ndeps step (A53).
.SFA_HESS_TIE_TOL <- 1e-8

## The Hessian for a point substituted for optim()'s, from the estimator
## optim(hessian = TRUE) itself uses, so standard errors do not depend on which
## stage won; at a tie, optim()'s own Hessian is carried across. A53, see
## notes/code_history/opts.md. `control` must match the final optim() call's.
.sfa_substitute_hessian <- function(fn, par, ref = NULL, control = list()) {
  rp <- ref$par
  rh <- ref$hessian
  np <- length(par)
  if (!is.null(rp) && length(rp) == np && all(is.finite(rp)) &&
    is.matrix(rh) && all(dim(rh) == np) && all(is.finite(rh)) &&
    all(abs(par - rp) <= .SFA_HESS_TIE_TOL * pmax(1, abs(rp)))) {
    return(rh)
  }
  h <- tryCatch(stats::optimHess(par, fn, control = control), error = function(e) NULL)
  if (is.null(h) || !all(is.finite(h))) NULL else h
}

## Stage 3 is a refinement pass over a point stages 1-2 already produced, so a
## failure here should cost the polish, not the fit. See
## notes/code_history/opts.md for the two failures this guards.
##
## NOT here, deliberately: a guard that falls back to the stage-2 point whenever
## optim() ends at a higher objective than it started from. It looks obviously
## right -- stages 1 and 2 refuse to move backwards, so why not stage 3 -- and it
## is wrong, because "the point it started from" is not always trustworthy. On
## sfma()'s own test data the NR likelihood has a spike: the stage-2 point
## evaluates to -3.6e17, optim() correctly escapes it to 541.8, and such a guard
## drags the fit back onto the spike and hands sfma() a log-likelihood of
## +3.6e17. L-BFGS-B is a descent method; when it ends above its own start, that
## is evidence about the START, not about the search. See
## test-ttsfm-numerics.R, which pins this.
## `gr` and `hess_fn` are optional exact derivatives (psfm()'s derivatives =
## "autograd", gap N3): gr goes to optim(), and hess_fn replaces the
## finite-difference Hessian, at the returned point and at an A42 substitute.
opt.optim <- function(fn, start_v, lower.optim, upper.optim, maxit.optim, opt.TF, method, optHessian, trace, verbose = verbose,
                      gr = NULL, hess_fn = NULL) {
  ## This function always hands optim() bounds, so the method must be one that
  ## honours them (gap A57).
  .check_optim_method(method)
  start_feval <- fn(start_v)
  opt <- NULL
  if (isTRUE(opt.TF == TRUE)) {
    np <- length(start_v)

    run_optim <- function(use_hessian) {
      tryCatch(
        optim(
          par = start_v,
          fn = fn,
          gr = gr,
          lower = lower.optim,
          upper = upper.optim,
          hessian = use_hessian,
          method = method,
          control = list(
            maxit = maxit.optim,
            REPORT = base::ceiling(maxit.optim / 10),
            trace = if (verbose) {
              1
            } else {
              0
            }
          )
        ),
        error = function(e) NULL
      )
    }

    opt <- run_optim(isTRUE(optHessian) && is.null(hess_fn))
    if (!is.null(opt) && isTRUE(optHessian) && !is.null(hess_fn) &&
      is.numeric(opt$value) && is.finite(opt$value)) {
      opt$hessian <- tryCatch(hess_fn(opt$par), error = function(e) NULL)
    }

    ## optim() threw, or stopped at a point the objective cannot evaluate.
    if (!is.null(opt) && (!is.numeric(opt$value) || !is.finite(opt$value))) {
      opt <- NULL
    }

    ## The point is usable but the numerical Hessian is not. Keep the point and
    ## hand downstream an all-NA Hessian, which every SE path already degrades
    ## to NA on.
    if (!is.null(opt) && isTRUE(optHessian) &&
      (is.null(opt$hessian) || !all(is.finite(opt$hessian)))) {
      opt$hessian <- matrix(NA_real_, np, np)
      opt$convergence <- .SFA_OPTIM_SKIPPED
      opt$message <- "final optim() stage returned a non-finite Hessian; standard errors are unavailable"
    }

    ## Nothing usable came back: return an optim-shaped result at the stage-2
    ## point so callers that index opt$par / opt$hessian keep working.
    if (is.null(opt)) {
      opt <- list(
        par = start_v,
        value = start_feval,
        counts = c(`function` = NA_integer_, gradient = NA_integer_),
        convergence = .SFA_OPTIM_SKIPPED,
        message = "final optim() stage skipped: non-finite objective or optim() error; the stage-2 result is returned",
        hessian = if (isTRUE(optHessian)) matrix(NA_real_, np, np) else NULL
      )
    }

    ## convergence is visible in print()/summary(), but coef() and friends do
    ## not go through those, so say it once out loud as well.
    if (identical(opt$convergence, .SFA_OPTIM_SKIPPED)) {
      warning(
        "The final optim() stage did not produce a usable result (", opt$message,
        "). Estimates come from the preceding stage and standard errors are NA. ",
        "Try optHessian = FALSE, vcov(type = \"bhhh\"), or different starting values.",
        call. = FALSE
      )
    }

    if (isTRUE(start_feval > opt$value)) {
      start_v <- opt$par
      start_feval <- fn(start_v)
    }

    ## A42. start_v/start_feval now hold the better of the two, but callers
    ## report `opt`, so when optim() ends ABOVE the point it was handed the
    ## better answer is computed and then discarded. Taking it unconditionally
    ## is the known-wrong fix described above, so it is taken only when the
    ## point looks like an optimum rather than a boundary collapse. See
    ## notes/code_history/opts.md.
    if (!is.null(opt) && is.finite(start_feval) &&
      is.numeric(opt$value) && is.finite(opt$value) &&
      opt$value > start_feval + 1e-10 &&
      .sfa_point_credible(fn, start_v, start_feval, lower.optim, upper.optim)) {
      h <- if (isTRUE(optHessian) && !is.null(hess_fn)) {
        tryCatch(hess_fn(start_v), error = function(e) NULL)
      } else if (isTRUE(optHessian)) {
        .sfa_substitute_hessian(fn, start_v, ref = opt)
      } else {
        NULL
      }
      ## The reported point must carry its own Hessian, or the standard errors
      ## would describe a point that is not the one being reported.
      if (!isTRUE(optHessian) || (!is.null(h) && all(is.finite(h)))) {
        opt$par <- start_v
        opt$value <- start_feval
        if (isTRUE(optHessian)) opt$hessian <- h
      }
    }
  }

  results <- list(start_v, start_feval, opt)
  names(results) <- c("start_v", "start_feval", "opt")
  return(results)
}

opt.psoptim <- function(fn, start_v, lower.psoptim, upper.psoptim = NA, maxit.psoptim = maxit.psoptim,
                        psopt.TF, rand.order = TRUE, verbose = verbose, rand.psoptim = rand.psoptim) {
  start_feval <- fn(start_v)
  opt00 <- NULL
  report <- base::ifelse(verbose, base::ceiling(maxit.psoptim / 10), 0)
  if (isTRUE(psopt.TF == TRUE)) {
    if (!is.null(rand.psoptim)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand.psoptim)
    }

    opt00 <- psoptim(
      par = start_v,
      fn = fn,
      lower = lower.psoptim,
      upper = upper.psoptim,
      control = list(
        trace = if (verbose) {
          1
        } else {
          0
        },
        REPORT = report,
        trace.stats = if (verbose) {
          TRUE
        } else {
          FALSE
        },
        maxit = maxit.psoptim,
        rand.order = rand.order
      )
    )

    if (isTRUE(start_feval > opt00$value)) {
      start_v <- opt00$par
      start_feval <- fn(start_v)
    }
  }

  results <- list(start_v, start_feval, opt00)
  names(results) <- c("start_v", "start_feval", "opt00")
  return(results)
}

## Multistart Newton with exact derivatives, the alternative to the bobyqa ->
## psoptim -> optim pipeline (psfm(optimizer = "newton"), gap N3). nlminb()
## from every start with the exact gradient and Hessian; the lowest objective
## wins. Several starts because one is not enough: on GTRE the likelihood has
## a second mode where sigma_h collapses to 0, and a single Newton run from a
## naive start found it (notes/code_history/opts.md). Returns an optim-shaped
## `opt`, or NULL if no start produced a finite objective.
.sfa_newton_multistart <- function(fn, gr, he, starts, lower, maxit, optHessian) {
  runs <- lapply(starts, function(s0) {
    tryCatch(
      suppressWarnings(stats::nlminb(pmax(s0, lower), function(x) fn(x), gr, he,
        lower = lower, control = list(iter.max = maxit, eval.max = 2L * maxit)
      )),
      error = function(e) NULL
    )
  })
  obj <- vapply(runs, function(r) {
    if (is.null(r) || !is.numeric(r$objective) || !is.finite(r$objective)) NA_real_ else r$objective
  }, numeric(1))
  tab <- data.frame(
    start = seq_along(starts), objective = obj,
    convergence = vapply(runs, function(r) if (is.null(r)) NA_integer_ else as.integer(r$convergence), integer(1)),
    iterations = vapply(runs, function(r) if (is.null(r)) NA_integer_ else as.integer(r$iterations), integer(1))
  )
  if (all(is.na(obj))) return(NULL)
  b <- which.min(obj)
  r <- runs[[b]]
  np <- length(r$par)
  h <- NULL
  if (isTRUE(optHessian)) {
    h <- tryCatch(he(r$par), error = function(e) NULL)
    if (is.null(h) || !all(is.finite(h))) h <- matrix(NA_real_, np, np)
  }
  list(
    opt = list(
      par = r$par, value = r$objective,
      counts = c(`function` = unname(r$evaluations[1]), gradient = unname(r$evaluations[2])),
      convergence = as.integer(r$convergence), message = r$message, hessian = h
    ),
    starts = tab, chosen = b, start_points = do.call(rbind, starts)
  )
}

## Starting points for .sfa_newton_multistart() on GTRE/TRE: the package's own
## random-effects start first, then the persistent scale spread over its
## split (GTRE: sigma_r vs sigma_h at three angles and two sizes) or over its
## size (TRE). Layout (lambda, sigma, sigma_r[, sigma_h], beta).
.sfa_newton_starts_panel <- function(start_v, model_name) {
  sv <- start_v
  sig <- if (is.finite(sv[2]) && sv[2] > 0) sv[2] else 1
  if (identical(model_name, "GTRE")) {
    s <- sqrt(sv[3]^2 + sv[4]^2)
    if (!is.finite(s) || s < 0.05 * sig) s <- 0.25 * sig
    grid <- list(c(1, 1), c(1, 2), c(1, 3), c(0.5, 2), c(2, 2))
    alt <- lapply(grid, function(a) {
      x <- sv
      x[3] <- a[1] * s * cos(a[2] * pi / 8)
      x[4] <- a[1] * s * sin(a[2] * pi / 8)
      x
    })
  } else {
    s <- if (is.finite(sv[3]) && sv[3] > 0.05 * sig) sv[3] else 0.25 * sig
    alt <- lapply(c(0.5, 2), function(m) {
      x <- sv
      x[3] <- m * s
      x
    })
  }
  c(list(sv), alt)
}
