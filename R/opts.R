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

opt.optim <- function(fn, start_v, lower.optim, upper.optim, maxit.optim, opt.TF, method, optHessian, trace, verbose = verbose) {
  start_feval <- fn(start_v)
  opt <- NULL
  if (isTRUE(opt.TF == TRUE)) {
    opt <- optim(
      par = start_v,
      fn = fn,
      lower = lower.optim,
      upper = upper.optim,
      hessian = optHessian,
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
    )

    if (isTRUE(start_feval > opt$value)) {
      start_v <- opt$par
      start_feval <- fn(start_v)
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
