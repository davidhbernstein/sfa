lcsfm <- function(formula,
                  model_name = c("LCM", "LCM_Z", "LCM_CN"),
                  data,
                  n_class = 2,
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
                  penalty_c = 0,
                  upper = NA,
                  Method = "L-BFGS-B",
                  verbose = FALSE,
                  keep_objective = FALSE,
                  rand.psoptim = NULL) {
  ## call/model_name resolution ahead of .check_model_formula_pipes() -- see
  ## sfm.R's identical fix for why.
  call <- match.call()
  model_name <- .match_model_name(model_name, eval(formals()$model_name))

  if (length(n_class) != 1L || !is.finite(n_class) ||
    n_class != as.integer(n_class) || n_class < 2L) {
    stop("lcsfm(): `n_class` must be a single integer >= 2; got ",
      paste(deparse(n_class), collapse = " "), ".",
      call. = FALSE
    )
  }
  n_class <- as.integer(n_class)

  .validate_sfa_call(formula, data, "lcsfm",
    maxit = list(
      maxit.bobyqa = maxit.bobyqa, maxit.psoptim = maxit.psoptim,
      maxit.optim = maxit.optim
    ),
    flags = list(optHessian = optHessian, PSopt = PSopt, inefdec = inefdec)
  )

  .check_model_formula_pipes(formula, model_name)

  DR1 <- data_proc(formula, data, model_name, individual = NULL, inefdec)

  data_orig <- DR1$data_orig
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

  Start_Cs <- start_cs(formula_x, data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val, n_z_vars, z_vars, n_class = n_class)

  lower_bob <- Start_Cs$lower_bob
  out <- Start_Cs$out
  start_v <- Start_Cs$start_v

  DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num = NA, individual = NA, N, model_name, rand.gtre = NULL)

  data <- DR2$data
  Y <- DR2$Y
  data_i_vars <- DR2$data_i_vars

  ## Latent class frontier (Greene 2005; Orea and Kumbhakar 2004; Caudill 2003): J
  ## technologies, class unobserved, log L_i = log sum_j P(j | q_i) f_j(eps_ij). "LCM"
  ## fixes the class probabilities, "LCM_Z" models them. See notes/code_history/lcsfm.md.
  if (length(penalty_c) != 1L || !is.numeric(penalty_c) || !is.finite(penalty_c) ||
    penalty_c < 0) {
    stop("lcsfm(): `penalty_c` must be a single finite number >= 0.",
      call. = FALSE
    )
  }
  if (penalty_c > 0 && identical(model_name, "LCM_Z")) {
    stop("lcsfm(): `penalty_c` applies to the unconditional class ",
      "probabilities, which \"LCM_Z\" makes depend on covariates. The ",
      "modified-likelihood result of Chen et al. (2001) is stated for a ",
      "scalar mixing proportion, so it is offered for \"LCM\" only.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------------
  ## LCM_CN, the contaminated normal frontier: only the noise scale differs across
  ## components, so f_eps is a closed-form mixture of NHN densities sharing sigma_u --
  ## the specification lcsfm_homogeneity()'s chi^2_{0:1} null is established for.
  ## ---------------------------------------------------------------------
  if (identical(model_name, "LCM_CN")) {
    J <- n_class
    Xm <- as.matrix(data_i_vars)
    n_b <- n_x_vars
    ## Layout: [sigv_1..sigv_J, sigu, beta(1..n_b), logit(1..J-1)].
    .lcm_pos <- c(rep(TRUE, J + 1L), rep(FALSE, n_b + (J - 1L)))

    .cn_logf <- function(x) {
      sigv <- pmax(abs(x[seq_len(J)]), .SFA_CONSTANTS$MIN_POSITIVE)
      sigu <- max(abs(x[J + 1L]), .SFA_CONSTANTS$MIN_POSITIVE)
      b <- x[(J + 2L):(J + 1L + n_b)]
      eps <- inefdec_n * (Y - Xm %*% b)
      lf <- matrix(0, nrow = length(Y), ncol = J)
      for (j in seq_len(J)) {
        sig <- sqrt(sigv[j]^2 + sigu^2)
        lam <- sigu / sigv[j]
        z2 <- pmin(pmax(-eps * lam / sig,
          .SFA_CONSTANTS$CLIP_Z1_LOWER
        ), .SFA_CONSTANTS$CLIP_Z1_UPPER)
        lf[, j] <- log(2) - log(sig) +
          stats::dnorm(eps / sig, log = TRUE) +
          stats::pnorm(z2, log.p = TRUE)
      }
      lf
    }

    ## Constant class probabilities: the contamination share is a scalar, which
    ## is the whole point of the restriction.
    .cn_logpi <- function(x) {
      d0 <- J + 1L + n_b
      eta <- matrix(0, nrow = length(Y), ncol = J)
      if (J > 1L) {
        for (j in seq_len(J - 1L)) {
          eta[, j] <- pmin(pmax(x[d0 + j],
            -.SFA_CONSTANTS$EXP_CLIP_UPPER
          ), .SFA_CONSTANTS$EXP_CLIP_UPPER)
        }
      }
      eta - .log_row_sum_exp(eta)
    }

    .cn_penalty <- function(x) {
      if (!isTRUE(penalty_c > 0)) return(0)
      penalty_c * (J * log(J) + sum(.cn_logpi(x)[1L, ]))
    }

    like.fn <- function(x, per_obs = FALSE) {
      like <- .log_row_sum_exp(.cn_logpi(x) + .cn_logf(x))
      like[!is.finite(like)] <- -sqrt(.Machine$double.xmax / length(like))
      if (isTRUE(per_obs)) return(like)
      -(sum(like[is.finite(like)]) + .cn_penalty(x))
    }

    Start.Time <- start.time()
    Opt.Bobyqa <- opt.bobyqa(
      fn = like.fn, start_v = start_v, lower.bobyqa = lower_bob,
      maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    bob1 <- Opt.Bobyqa$bob1

    attr(start_v, "lcm_pos") <- .lcm_pos
    Lower.Start <- lower.start(start_v, "LCM", differ = 1)
    Opt.Psoptim <- opt.psoptim(
      fn = like.fn, start_v, lower.psoptim = Lower.Start$lower1,
      rand.psoptim = rand.psoptim, upper.psoptim = Lower.Start$upper1,
      maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    opt00 <- Opt.Psoptim$opt00

    attr(start_v, "lcm_pos") <- .lcm_pos
    Lower.Start <- lower.start(start_v, "LCM", differ = 0.5)
    Opt.Optim <- opt.optim(
      fn = like.fn, start_v = start_v, lower.optim = Lower.Start$lower1,
      upper.optim = Lower.Start$upper1_open, maxit.optim = maxit.optim,
      opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    opt <- Opt.Optim$opt
    End.Time <- end.time(Start.Time)

    ## The scales enter through abs(), so their sign is not identified and the
    ## optimizer may return either. Report the magnitude.
    .scale_at <- which(.lcm_pos)
    start_v[.scale_at] <- abs(start_v[.scale_at])
    st_err <- if (is.null(opt$hessian) || all(!is.finite(opt$hessian)) ||
      isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA_real_, length(start_v))
    } else {
      suppressWarnings(sqrt(diag(solve(opt$hessian))))
    }
    out[1, ] <- start_v
    out[2, ] <- st_err
    out[3, ] <- out[1, ] / out[2, ]
    rownames(out) <- c("par", "st_err", "t-val")

    lpi <- .cn_logpi(start_v)
    class_prob <- colMeans(exp(lpi))
    names(class_prob) <- paste0("class", seq_len(J))
    post <- exp(lpi + .cn_logf(start_v) -
      .log_row_sum_exp(lpi + .cn_logf(start_v)))
    colnames(post) <- paste0("class", seq_len(J))
    .pen <- .cn_penalty(start_v)

    results <- list(
      t(out), c(opt), End.Time, start_v, model_name, formula, post, J,
      ## MINUS .pen, not plus -- see the identical correction in the
      ## LCM/LCM_Z block below. like.fn returns -(L + P), so -opt$value is
      ## L + P and the plain log-likelihood is -opt$value - P.
      class_prob, .pen, -opt$value - .pen, penalty_c,
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "model_name", "formula",
      "post.prob", "n_class", "class_prob",
      "penalty", "logLik_unpenalised", "penalty_c",
      "coefficients", "std.errors", "t.values", "call"
    )
    ## Rows actually used, not rows supplied: bread() scales by this.
    results$nobs <- length(as.numeric(Y))
    ## Optionally retain the objective, so estfun()/vcov(type = "bhhh") can
    ## difference it into a per-observation score matrix after the fact.
    ## NOTE for LCM_CN: like.fn's per_obs branch returns the UNPENALISED
    ## per-observation contributions -- a penalty on the class probabilities
    ## is not a per-observation quantity. With the default penalty_c = 0 the
    ## two coincide exactly; with penalty_c > 0 the scores are those of the
    ## unpenalised likelihood, so the OPG is not the penalised estimator's.
    if (isTRUE(keep_objective)) results$objective <- like.fn
    return(results)
  }

  if (model_name %in% c("LCM", "LCM_Z")) {
    J <- n_class
    blk <- 2L + n_x_vars ## sigv, sigu, beta for one class
    Xm <- as.matrix(data_i_vars)
    Qm <- if (model_name == "LCM") {
      matrix(1, nrow = length(Y), ncol = 1L)
    } else {
      as.matrix(data_z)
    }
    n_q <- ncol(Qm)
    ## Which entries of the parameter vector are scales, and so bounded below
    ## by zero. lower.start() reads this off the attribute rather than by index,
    ## because the layout depends on J and on the number of regressors.
    .lcm_pos <- c(
      rep(c(TRUE, TRUE, rep(FALSE, n_x_vars)), J),
      rep(FALSE, (J - 1L) * n_q)
    )

    ## Class-wise log densities: an N x J matrix of log f_j(eps_ij).
    .lcm_logf <- function(x) {
      lf <- matrix(0, nrow = length(Y), ncol = J)
      for (j in seq_len(J)) {
        i0 <- (j - 1L) * blk
        sigv <- abs(x[i0 + 1L])
        sigu <- abs(x[i0 + 2L])
        b_j <- x[(i0 + 3L):(i0 + blk)]
        sigv <- max(sigv, .SFA_CONSTANTS$MIN_POSITIVE)
        sigu <- max(sigu, .SFA_CONSTANTS$MIN_POSITIVE)
        sig <- sqrt(sigv^2 + sigu^2)
        lam <- sigu / sigv
        eps <- inefdec_n * (Y - Xm %*% b_j)
        z2 <- pmin(pmax(-eps * lam / sig,
          .SFA_CONSTANTS$CLIP_Z1_LOWER
        ), .SFA_CONSTANTS$CLIP_Z1_UPPER)
        ## log.p / log = TRUE rather than log(pnorm(.)): the far tail of the
        ## half-normal factor underflows to zero on the natural scale for
        ## exactly the observations a badly separated class is trying to
        ## explain, and log(0) would discard them.
        lf[, j] <- log(2) - log(sig) +
          stats::dnorm(eps / sig, log = TRUE) +
          stats::pnorm(z2, log.p = TRUE)
      }
      lf
    }

    ## Multinomial-logit log class probabilities, class J the reference.
    .lcm_logpi <- function(x) {
      eta <- matrix(0, nrow = length(Y), ncol = J)
      if (J > 1L) {
        d0 <- J * blk
        for (j in seq_len(J - 1L)) {
          d_j <- x[(d0 + (j - 1L) * n_q + 1L):(d0 + j * n_q)]
          eta[, j] <- pmin(pmax(Qm %*% d_j,
            -.SFA_CONSTANTS$EXP_CLIP_UPPER
          ), .SFA_CONSTANTS$EXP_CLIP_UPPER)
        }
      }
      eta - .log_row_sum_exp(eta)
    }

    ## Chen et al. (2001)'s penalty c log(J^J prod p_j), for L3's modified LR test: keeps
    ## every class probability off 0 and 1, and is exactly 0 at equal probabilities. "LCM" only.
    .lcm_penalty <- function(x) {
      if (!isTRUE(penalty_c > 0)) return(0)
      lp <- .lcm_logpi(x)[1L, ]
      penalty_c * (J * log(J) + sum(lp))
    }

    like.fn <- function(x, per_obs = FALSE) {
      like <- .log_row_sum_exp(.lcm_logpi(x) + .lcm_logf(x))
      like[!is.finite(like)] <- -sqrt(.Machine$double.xmax / length(like))
      if (isTRUE(per_obs)) return(like)
      -(sum(like[is.finite(like)]) + .lcm_penalty(x))
    }

    Start.Time <- start.time()

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.fn, start_v = start_v, lower.bobyqa = lower_bob,
      maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    bob1 <- Opt.Bobyqa$bob1

    ## lower.start() needs the block structure back; the optimizer stages return
    ## a bare numeric and drop attributes, so re-attach it at each call.
    attr(start_v, "lcm_pos") <- .lcm_pos
    Lower.Start <- lower.start(start_v, model_name, differ = 1)
    Opt.Psoptim <- opt.psoptim(
      fn = like.fn, start_v, lower.psoptim = Lower.Start$lower1,
      rand.psoptim = rand.psoptim, upper.psoptim = Lower.Start$upper1,
      maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    opt00 <- Opt.Psoptim$opt00

    attr(start_v, "lcm_pos") <- .lcm_pos
    Lower.Start <- lower.start(start_v, model_name, differ = 0.5)
    Opt.Optim <- opt.optim(
      fn = like.fn, start_v = start_v, lower.optim = Lower.Start$lower1,
      upper.optim = Lower.Start$upper1_open, maxit.optim = maxit.optim,
      opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- bob1
      st_err <- rep(NA, length(opt$par))
    }
    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
      st_err <- rep(NA, length(opt$par))
    }
    if (optHessian == TRUE) {
      st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
        rep(NA, length(opt$par))
      } else {
        suppressWarnings(sqrt(diag(solve(opt$hessian))))
      }
    }
    ## The two scales per class enter the likelihood through abs(), so their
    ## sign is not identified and the optimizer may return either. Report the
    ## magnitude, which is what the model is about.
    .scale_at <- which(.lcm_pos)
    opt$par[.scale_at] <- abs(opt$par[.scale_at])

    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## Posterior class probabilities P(j | i), which weight the class-conditional JLMS
    ## scores: the firm's frontier is uncertain, so no single modal-class score is reported.
    lf <- .lcm_logf(opt$par)
    lpi <- .lcm_logpi(opt$par)
    log_f <- .log_row_sum_exp(lpi + lf)
    post.prob <- exp((lpi + lf) - log_f)
    colnames(post.prob) <- paste0("class", seq_len(J))

    jlms_class <- matrix(0, nrow = length(Y), ncol = J)
    for (j in seq_len(J)) {
      i0 <- (j - 1L) * blk
      sigv <- max(abs(opt$par[i0 + 1L]), .SFA_CONSTANTS$MIN_POSITIVE)
      sigu <- max(abs(opt$par[i0 + 2L]), .SFA_CONSTANTS$MIN_POSITIVE)
      b_j <- opt$par[(i0 + 3L):(i0 + blk)]
      eps <- inefdec_n * (Y - Xm %*% b_j)
      sigsq <- sigv^2 + sigu^2
      mustar <- -eps * sigu^2 / sigsq
      sigstar <- sqrt(sigu^2 * sigv^2 / sigsq)
      zz <- mustar / sigstar
      jlms_class[, j] <- .jlms_u(mustar, sigstar)
    }
    colnames(jlms_class) <- paste0("class", seq_len(J))
    jlms <- rowSums(post.prob * jlms_class)
    class_assign <- max.col(post.prob, ties.method = "first")
    ## Prior (unconditional) class shares, averaged over the sample -- the
    ## number a reader wants when asking how big each technology group is.
    class_prob <- colMeans(exp(lpi))
    names(class_prob) <- paste0("class", seq_len(J))

    ## With penalty_c > 0 the optimizer maximised the MODIFIED likelihood, so
    ## opt$value is the penalised objective. logLik() must not report that as
    ## the log-likelihood, so the penalty and the plain log-likelihood are
    ## carried separately; lcsfm_homogeneity() needs both, and .sfa_penalty is
    ## exactly 0 (not merely small) whenever penalty_c is 0.
    ## At opt$par, which is what out[1, ] reports and what post.prob is built
    ## from -- not start_v, which the optimizer stages overwrite.
    .pen <- .lcm_penalty(opt$par)
    results <- list(
      t(out), c(opt), End.Time, start_v, model_name, formula, jlms, post.prob,
      jlms_class, class_assign, class_prob, J,
      ## MINUS .pen, not plus. like.fn returns -(L + P), so -opt$value is
      ## L + P and the plain log-likelihood is -opt$value - P. Adding it gave
      ## L + 2P, double-counting the penalty downward (P <= 0 always, by
      ## Jensen), so the stored "unpenalised" figure was further from the
      ## truth than the penalised objective it was meant to correct.
      .pen, -opt$value - .pen, penalty_c,
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "model_name", "formula", "jlms",
      "post.prob", "jlms_class", "class", "class_prob", "n_class",
      "penalty", "logLik_unpenalised", "penalty_c",
      "coefficients", "std.errors", "t.values", "call"
    )
    ## Rows actually used, not rows supplied: bread() scales by this.
    results$nobs <- length(as.numeric(Y))
    ## Optionally retain the objective, so estfun()/vcov(type = "bhhh") can
    ## difference it into a per-observation score matrix after the fact.
    ## NOTE for LCM_CN: like.fn's per_obs branch returns the UNPENALISED
    ## per-observation contributions -- a penalty on the class probabilities
    ## is not a per-observation quantity. With the default penalty_c = 0 the
    ## two coincide exactly; with penalty_c > 0 the scores are those of the
    ## unpenalised likelihood, so the OPG is not the penalised estimator's.
    if (isTRUE(keep_objective)) results$objective <- like.fn
    return(results)
  }


  stop(paste0(
    "model_name '", model_name, "' is a recognized choice for lcsfm() but has ",
    "no implementation branch. Valid choices are: \"LCM\", \"LCM_Z\"."
  ), call. = FALSE)
}
