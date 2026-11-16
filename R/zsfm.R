zsfm <- function(formula,
                 model_name = c("ZISF", "ZISF_Z"),
                 data,
                 maxit.bobyqa = 10000,
                 maxit.psoptim = 1000,
                 maxit.optim = 1000,
                 REPORT = 1,
                 trace = 0,
                 pgtol = 0,
                 start_val = FALSE,
                 PSopt = FALSE,
                 optHessian = TRUE,
                 inefdec = TRUE,
                 upper = NA,
                 Method = "L-BFGS-B",
                 logit = TRUE,
                 verbose = FALSE,
                 keep_objective = FALSE,
                 rand.psoptim = NULL,
                 ghet = NULL,
                 vhet = NULL) {
  ## call/model_name resolution moved ahead of .check_model_formula_pipes() --
  ## see sfm.R's identical fix for why.
  call <- match.call()
  data_supplied <- if (!missing(data)) data
  model_name <- .match_model_name(model_name, eval(formals()$model_name))

  ## HETEROSKEDASTIC NOISE, as a named argument rather than a further pipe
  ## segment. In zsfm() the reason is sharper than elsewhere: the second pipe
  ## part of ZISF_Z is the MIXING PROBABILITY's linear predictor, not a
  ## variance at all, so a reader cannot be asked to tell sigma_v from the
  ## regime link by position.
  ##
  ## Both branches of the mixture are already written elementwise in eps --
  ## the efficient regime is N(0, sigma_v^2) and the inefficient one is the
  ## ordinary normal/half-normal composite -- so sigma_v varying by
  ## observation changes what it IS and nothing about how either density is
  ## formed. The mixing probability is untouched: `ghet` and `vhet` are
  ## different blocks and may be given together.
  ##
  ## WITH `vhet` NULL NOTHING BELOW CHANGES.
  if (!is.null(vhet) && !inherits(vhet, "formula")) {
    stop("zsfm(): `vhet` must be a one-sided formula, e.g. ~ z1 + z2.",
      call. = FALSE
    )
  }

  .validate_sfa_call(formula, data, "zsfm",
    maxit = list(
      maxit.bobyqa = maxit.bobyqa, maxit.psoptim = maxit.psoptim,
      maxit.optim = maxit.optim
    ),
    flags = list(optHessian = optHessian, PSopt = PSopt, inefdec = inefdec)
  )

  ## A determinant block may be written as a pipe segment or as a named
  ## argument, and the named one is spliced into its pipe position here, so
  ## everything below this line sees one canonical formula and no likelihood,
  ## data_proc() tier or bounds section had to change.
  ##
  ## PIPES ARE UNCHANGED. With no named argument given this returns `formula`
  ## itself -- the identical object -- so an existing pipe call does not enter
  ## the new path at all and cannot be affected by it.
  formula <- .sfa_splice_het_formula(formula, model_name, list(ghet = ghet),
                                     "zsfm")
  .check_model_formula_pipes(formula, model_name)
  .check_optim_method(Method)

  ## Complete cases over `vhet` BEFORE data_proc(), so its row bookkeeping and
  ## the Zv design cannot fall out of step -- same ordering as sfm().
  if (!is.null(vhet)) data <- .het_prefilter(data, list(vhet))

  DR1 <- data_proc(formula, data, model_name, individual = NULL, inefdec)

  data_orig <- DR1$data_orig
  data <- data_orig ## carries a transformed response as its own column
  form_parts <- DR1$form_parts
  formula_x <- DR1$formula_x
  y_var <- DR1$y_var
  model_name <- DR1$model_name
  data_x <- DR1$data_x
  intercept <- DR1$intercept
  inefdec_n <- DR1$inefdec_n
  inefdec_TF <- DR1$inefdec_TF
  x_vars_vec <- DR1$x_vars_vec
  n_x_vars <- DR1$n_x_vars
  x_vars <- DR1$x_vars
  x_x_vec <- DR1$x_x_vec
  fancy_vars <- DR1$fancy_vars
  fancy_vars_z <- DR1$fancy_vars_z
  n_z_vars <- DR1$n_z_vars
  N <- DR1$N
  data_z <- DR1$data_z
  if (length(unlist(form_parts)) > 3) {
    formula_z <- DR1$formula_z
    intercept_z <- DR1$intercept_z
    n_z_vars <- DR1$n_z_vars
    z_vars <- DR1$z_vars
    z_vars_vec <- DR1$z_vars_vec
    z_z_vec <- DR1$z_z_vec
  }


  Start_Cs <- start_cs(formula_x, data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val, n_z_vars, z_vars)

  beta_0 <- Start_Cs$beta_0
  beta_0_st <- Start_Cs$beta_0_st
  beta_hat <- Start_Cs$beta_hat
  epsilon_hat <- Start_Cs$epsilon_hat
  lambda <- Start_Cs$lambda
  lower_bob <- Start_Cs$lower_bob
  mu <- Start_Cs$mu
  out <- Start_Cs$out
  plm_lm <- Start_Cs$plm_lm
  sigma <- Start_Cs$sigma
  sigma_u <- Start_Cs$sigma_u
  sigma_v <- Start_Cs$sigma_v
  start_v <- Start_Cs$start_v
  start_v_ne <- Start_Cs$start_v_ne
  start_v_ng <- Start_Cs$start_v_ng
  start_v_nhn <- Start_Cs$start_v_nhn
  start_v_nnak <- Start_Cs$start_v_nnak
  start_v_nr <- Start_Cs$start_v_nr
  start_v_ntn <- Start_Cs$start_v_ntn
  start_v_t <- Start_Cs$start_v_t

  DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num = NA, individual = NA, N, model_name, rand.gtre = NULL)

  data <- DR2$data
  Y <- DR2$Y
  data_i_vars <- DR2$data_i_vars

  ## The sigma_v design. Intercept-only when `vhet` is NULL, which is the
  ## homoskedastic case. Built after data_proc2(), the step that settles which
  ## rows survive.
  Zv <- .het_design(vhet, data, "vhet", caller = "zsfm")
  n_zv_vars <- ncol(Zv)

  ## sigma_v and its square from the block.
  ##
  ## The NULL branch reproduces the ORIGINAL expressions, x[i]^2 then sqrt(),
  ## rather than the algebraically equal abs(x[i]). sqrt(x^2) and abs(x) are
  ## not the same floating-point computation, and an existing fit must not
  ## move by a single bit because this argument was added.
  ##
  ## A given `vhet` puts the block on the LOG-SD scale, exp(z'delta). zsfm()
  ## has no `z_link` argument -- its one determinant block is the mixing
  ## probability, which is a logit, not a variance -- so there is no
  ## package-level link convention here to honour, and the log-SD form is the
  ## one sfm()'s default and ttsfm() both use.
  .zi_sigv <- function(d) {
    if (is.null(vhet)) {
      sigvsq <- d[1L]^2
      list(sigv = sqrt(sigvsq), sigvsq = sigvsq)
    } else {
      sv <- exp(pmin(as.numeric(Zv %*% d), .SFA_CONSTANTS$EXP_CLIP_UPPER))
      list(sigv = sv, sigvsq = sv^2)
    }
  }

  ## Where the sigma_v block sits: ZISF's layout is (gamma, sigv, sigu, beta),
  ## ZISF_Z's is (sigv, sigu, beta, gamma-block). sigma_u follows the block.
  .i_sv <- if (identical(model_name, "ZISF")) 2L else 1L
  .i_su <- .i_sv + n_zv_vars

  ## Splice the block into the layout start_cs() built. Done HERE, not in
  ## start_cs(), because that helper is shared with sfm() and every other
  ## cross-sectional model; widening one slot there would reach models that
  ## have no `vhet`.
  if (!is.null(vhet)) {
    .nm0 <- colnames(out)
    .blk <- .het_block_start(Zv, start_v[.i_sv], "sd")
    .splice <- function(v, new) append(v[-.i_sv], new, after = .i_sv - 1L)
    start_v <- .splice(start_v, .blk)
    ## The homoskedastic slot was bounded positive because it held sigma_v
    ## directly. The block holds a log-scale linear predictor, which is
    ## unbounded by construction, so the bound is dropped rather than widened.
    lower_bob <- .splice(lower_bob, rep(-Inf, n_zv_vars))
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- .splice(.nm0, paste0("Zv.", colnames(Zv)))
    rownames(out) <- c("par", "st_err", "t-val")
  }

  if (model_name %in% c("ZISF", "ZISF_Z")) {
    ## `per_obs = TRUE` returns the vector of per-observation log-likelihood
    ## contributions instead of the negative sum, which is what
    ## estfun.sfareg() differences to build the score matrix behind
    ## vcov(type = "bhhh"). The optimizers call this with one argument.
    like.fn <- function(x, per_obs = FALSE) {
      ## Width of the sigma_v block; 1 unless `vhet` was given, in which case
      ## every index after it shifts by nzv - 1.
      nzv <- n_zv_vars
      if (model_name %in% c("ZISF")) {
        x_x_vec <- x[(nzv + 2L) + seq_len(n_x_vars)]
      }
      if (model_name %in% c("ZISF_Z")) {
        x_x_vec <- x[(nzv + 1L) + seq_len(n_x_vars)]
      }

      eps <- (inefdec_n * (Y - as.matrix(data_i_vars) %*% x_x_vec))

      if (model_name == "ZISF") {
        gamma <- x[1]
        ## exp(-|gamma|), so the likelihood is EXACTLY symmetric in gamma and
        ## the optimizer may legitimately return either sign.
        prob <- exp(-abs(gamma))
        .sv <- .zi_sigv(x[1L + seq_len(nzv)])
        sigv <- .sv$sigv
        sigvsq <- .sv$sigvsq
        sigusq <- x[1L + nzv + 1L]^2
        sigu <- sqrt(sigusq)

        lambda <- sigu / sigv
        sigsq <- sigvsq + sigusq
        sig <- sqrt(sigsq)

        f1 <- -0.5 * log(2 * pi * sigvsq) - (0.5 / sigvsq) * eps^2
        f2 <- log(2 / sig) + dnorm(-eps / sig, log = TRUE) + pnorm(-eps * lambda / sig, log.p = TRUE)
        ## Mixed on the LOG scale.
        like <- .log_add2(log(prob) + f1, log1p(-prob) + f2)
      }

      if (model_name == "ZISF_Z") {
        gamma <- x[(nzv + 1L + n_x_vars) + seq_len(n_z_vars)] ## lets put gammas last

        ## plogis(), not exp(eta)/(1+exp(eta)): the explicit form overflows to
        ## Inf/Inf = NaN once eta passes ~710.
        if (logit) {
          prob <- plogis(data_z %*% gamma)
        }
        ## NOTE: this is pnorm(eta)/(1+pnorm(eta)), which is bounded ABOVE BY
        ## 0.5 -- it is not the probit link.
        if (!logit) {
          prob <- pnorm(data_z %*% gamma) / (1 + pnorm(data_z %*% gamma))
        }

        .sv <- .zi_sigv(x[seq_len(nzv)])
        sigv <- .sv$sigv
        sigvsq <- .sv$sigvsq
        sigusq <- x[nzv + 1L]^2
        sigu <- sqrt(sigusq)

        lambda <- sigu / sigv
        sigsq <- sigvsq + sigusq
        sig <- sqrt(sigsq)

        f1 <- -0.5 * log(2 * pi * sigvsq) - (0.5 / sigvsq) * eps^2
        f2 <- log(2 / sig) + dnorm(-eps / sig, log = TRUE) + pnorm(-eps * lambda / sig, log.p = TRUE)
        like <- .log_add2(log(prob) + f1, log1p(-prob) + f2)
      }


      like[like == -Inf] <- -.SFA_CONSTANTS$NONFINITE_PENALTY / length(like)
      like[like == Inf] <- -.SFA_CONSTANTS$NONFINITE_PENALTY / length(like)
      like[is.nan(like)] <- -.SFA_CONSTANTS$NONFINITE_PENALTY / length(like)
      ## The barrier above leaves NA untouched; per_obs must return one row
      ## per observation, so anything still non-finite is floored here rather
      ## than dropped the way the summed branch drops it.
      if (isTRUE(per_obs)) {
        like[!is.finite(like)] <- -.SFA_CONSTANTS$NONFINITE_PENALTY / length(like)
        return(like)
      }

      return(-sum(like[is.finite(like)]))
    }

    Start.Time <- start.time()

    ## Multi-start over the mixing parameter, for the same reason NG has one.
    ##
    ## P(fully efficient) = exp(-gamma), so gamma -> Inf is the boundary where
    ## the mass point vanishes and the model collapses to an ordinary frontier.
    ## That boundary has its own local optimum, and the default start lands in
    ## its basin on a small fraction of samples: measured on the convergence
    ## design, 1 replication in 5000 returned gamma = 9.93 (p = 5e-5) where a
    ## start from anywhere in the interior reached gamma = 0.31 with a
    ## log-likelihood 188.7 points HIGHER. One such fit was enough to inflate
    ## that sample size's mean MSE 144-fold and flatten the whole log-MSE-on-
    ## log-n slope to zero, which is what made gamma read as non-convergent
    ## when four of the five sample sizes were textbook.
    ##
    ## Cheap: the extra starts differ only in gamma, and the objective is
    ## closed-form. The best attained objective wins, so this can only improve
    ## the fit -- a start that does worse is discarded.
    ## ZISF ONLY. ZISF_Z has no scalar gamma: its layout is
    ## (sigv, sigu, beta, z-coefficients) and the mixing is parameterised
    ## through the z block at the END, so start_v[1] there is sigma_v and
    ## perturbing it would move the noise scale rather than the mixture.
    if (identical(model_name, "ZISF") && isFALSE(is.numeric(start_val))) {
      .zi_cand <- lapply(c(0.1, 0.3, 0.7, 1.5), function(g) replace(start_v, 1L, g))
      .zi_cand <- c(list(start_v), .zi_cand)
      .zi_obj <- vapply(.zi_cand, function(z) {
        tryCatch({
          v <- like.fn(z)
          if (is.finite(v)) v else Inf
        }, error = function(e) Inf)
      }, numeric(1))
      if (any(is.finite(.zi_obj))) {
        .zi_fits <- Filter(Negate(is.null), lapply(order(.zi_obj)[seq_len(min(3L, sum(is.finite(.zi_obj))))],
          function(i) {
            tryCatch(suppressWarnings(stats::optim(.zi_cand[[i]], like.fn,
              method = "L-BFGS-B", lower = lower_bob,
              control = list(maxit = 200))), error = function(e) NULL)
          }))
        .zi_fits <- Filter(function(o) is.finite(o$value), .zi_fits)
        if (length(.zi_fits)) {
          start_v <- .zi_fits[[which.min(vapply(.zi_fits, function(o) o$value, numeric(1)))]]$par
        }
      }
    }

    ## ZISF_Z has the same problem through its regime link: a spurious basin
    ## about 40 log-likelihood units low (gap A31). Try several link intercepts
    ## with zero slopes, polish the best three, keep the best.
    if (identical(model_name, "ZISF_Z") && isFALSE(is.numeric(start_val))) {
      .gi <- (n_zv_vars + 1L + n_x_vars) + seq_len(n_z_vars)
      .zz_cand <- lapply(c(-2, -0.4, 0.4, 2), function(g0) {
        replace(start_v, .gi, c(g0, rep(0, length(.gi) - 1L)))
      })
      .zz_cand <- c(list(start_v), .zz_cand)
      .zz_obj <- vapply(.zz_cand, function(z) {
        tryCatch({
          v <- like.fn(z)
          if (is.finite(v)) v else Inf
        }, error = function(e) Inf)
      }, numeric(1))
      if (any(is.finite(.zz_obj))) {
        .zz_fits <- Filter(Negate(is.null), lapply(order(.zz_obj)[seq_len(min(3L, sum(is.finite(.zz_obj))))],
          function(i) {
            tryCatch(suppressWarnings(stats::optim(.zz_cand[[i]], like.fn,
              method = "L-BFGS-B", lower = lower_bob,
              control = list(maxit = 200))), error = function(e) NULL)
          }))
        .zz_fits <- Filter(function(o) is.finite(o$value), .zz_fits)
        if (length(.zz_fits)) {
          start_v <- .zz_fits[[which.min(vapply(.zz_fits, function(o) o$value, numeric(1)))]]$par
        }
      }
    }

    Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v, lower.bobyqa = lower_bob, maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose)
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    ## lower.start() lays its bound vectors out from the HOMOSKEDASTIC
    ## parameter order: it bounds two adjacent slots positive, because in that
    ## order they are sigma_v and sigma_u. With a sigma_v BLOCK in place of
    ## the first, both bounds land on the wrong slots -- the positive one on a
    ## LOG-scale predictor, which may legitimately be negative, and none at
    ## all on sigma_u.
    ##
    ## Pinning a log-scale slot at >= 1e-7 holds sigma_v at about 1 whatever
    ## the data says. That is not hypothetical: it is what the first run of
    ## this feature did, and it drove the mixing parameter to its boundary and
    ## left sigma_u with a standard error of 21.
    .ls_het <- function(differ) {
      lo <- start_v - differ
      lo[.i_su] <- .0000001
      up <- start_v + differ
      list(lower1 = lo, upper1 = up,
           upper1_open = replace(up, .i_su, Inf))
    }
    Lower.Start <- if (is.null(vhet)) {
      lower.start(start_v, model_name, differ = 1)
    } else {
      .ls_het(1)
    }

    Opt.Psoptim <- opt.psoptim(
      fn = like.fn, start_v, lower.psoptim = Lower.Start$lower1, rand.psoptim = rand.psoptim,
      upper.psoptim = Lower.Start$upper1, maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    Lower.Start <- if (is.null(vhet)) {
      lower.start(start_v, model_name, differ = 0.5)
    } else {
      .ls_het(0.5)
    }
    Opt.Optim <- opt.optim(
      fn = like.fn, start_v = start_v, lower.optim = Lower.Start$lower1,
      upper.optim = Lower.Start$upper1_open, maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)
    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- .as_optim(bob1)
      st_err <- rep(NA, length(opt$par))
    }

    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
      st_err <- rep(NA, length(opt$par))
    }

    if (optHessian == TRUE) {
      ## The colMeans test catches an all-zero Hessian only; a singular but
      ## non-zero one made solve() error and discard a finished fit (A65).
      st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
        rep(NA, length(opt$par))
      } else {
        tryCatch(suppressWarnings(sqrt(diag(solve(opt$hessian)))),
          error = function(e) rep(NA_real_, length(opt$par))
        )
      }
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## suppressWarnings() above discards R's own "NaNs produced", which was the
    ## only signal that a standard error had gone missing (A65).
    .sm <- if (optHessian == TRUE) {
      .se_missing_report(out["st_err", ], colnames(out), opt$hessian,
        paste0("zsfm(model_name = \"", model_name, "\"): ")
      )
    } else {
      NULL
    }
    if (!is.null(.sm)) warning(.sm$message, call. = FALSE)

    ## JLMS
    if (model_name %in% c("ZISF", "ZISF_Z")) {
      ## Branch on the MODEL, not on is.na(n_z_vars).
      nzv <- n_zv_vars
      if (model_name == "ZISF") {
        beta <- opt$par[-seq_len(nzv + 2L)]
        z <- 1
        gamma <- opt$par[1]
        ## exp(-|gamma|), matching the likelihood that was actually maximised.
        prob <- exp(-abs(gamma))

        .sv <- .zi_sigv(opt$par[1L + seq_len(nzv)])
        sigvsq <- .sv$sigvsq
        sigusq <- opt$par[1L + nzv + 1L]^2
      }

      if (model_name == "ZISF_Z") {
        beta <- opt$par[(nzv + 1L) + seq_len(n_x_vars)]
        gamma <- opt$par[(nzv + 1L + n_x_vars) + seq_len(n_z_vars)] ## lets put gammas last

        if (logit) {
          prob <- plogis(data_z %*% gamma)
        }
        if (!logit) {
          prob <- pnorm(data_z %*% gamma) / (1 + pnorm(data_z %*% gamma))
        }

        .sv <- .zi_sigv(opt$par[seq_len(nzv)])
        sigvsq <- .sv$sigvsq
        sigusq <- opt$par[nzv + 1L]^2
      }

      eps <- (inefdec_n * (Y - as.matrix(data_i_vars) %*% beta))
      sigv <- .sv$sigv
      sigu <- sqrt(sigusq)

      ## Reparametrize the log-likelihood function
      lambda <- sigu / sigv
      sigsq <- sigvsq + sigusq
      sig <- sqrt(sigsq)

      ## Now the likelihood function, on the log scale as above so that the
      ## posterior probability stays a ratio of quantities that cannot underflow.
      f1 <- -0.5 * log(2 * pi * sigvsq) - (0.5 / sigvsq) * eps^2
      f2 <- log(2 / sig) + dnorm(-eps / sig, log = TRUE) + pnorm(-eps * lambda / sig, log.p = TRUE)
      l_eff <- log(prob) + f1 ## efficient-regime contribution
      l_ineff <- log1p(-prob) + f2 ## inefficient-regime contribution
      log_f <- .log_add2(l_eff, l_ineff)

      post.prob <- exp(l_eff - log_f)
      mustar <- -eps * sigusq / sigsq
      sigstarsq <- sigusq * sigvsq / (sigusq + sigvsq)
      sigstar <- sqrt(sigstarsq)
      zz <- mustar / sigstar
      ## In logs: dnorm(zz) / pnorm(zz) is 0 / 0 below zz ~ -37 (gap A30).
      jlms <- .jlms_u(mustar, sigstar)
    }

    if (model_name %in% c("ZISF", "ZISF_Z")) {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, jlms, post.prob,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "jlms", "post.prob",
        "coefficients", "std.errors", "t.values", "call"
      )
      results$inefdec <- isTRUE(inefdec)
      ## Appended, not rebuilt: re-listing twelve values against twelve names
      ## to add one field is how a misalignment gets introduced.
      if (!is.null(.sm)) results$se_unresolved <- .sm
      ## Rows actually used, not rows supplied: bread() scales by this.
      results$nobs <- length(as.numeric(Y))
      results$anchor_resid <- .sfa_anchor_resid(data_i_vars, Y)
      ## Optionally retain the objective, so estfun()/vcov(type = "bhhh") and
      ## sfa_diagnostics() can difference it after the fact.
      if (isTRUE(keep_objective)) results$objective <- like.fn
      ## With `vhet` given there is no single sigma_v row in out[, "par"] --
      ## those rows are COEFFICIENTS -- so the fitted per-observation scale is
      ## stored rather than left to be rebuilt from the deltas.
      if (!is.null(vhet)) {
        results$vhet <- vhet
        results$sigma_v <- as.numeric(sigv)
      }
    }
    return(.sfa_keep_data(results, data_supplied))
  } else {
    stop(paste0(
      "model_name '", model_name, "' is a recognized choice for zsfm() but has no implementation branch. ",
      "Valid choices are: \"ZISF\", \"ZISF_Z\"."
    ), call. = FALSE)
  }
}
