ttsfm <- function(formula,
                  model_name = c("TTNE", "TTHN", "TTNLS"),
                  data,
                  z_link = c("sd", "var"),
                  maxit.bobyqa = 80000,
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
                  rand.psoptim = NULL) {
  ## call/model_name resolution moved ahead of .check_model_formula_pipes() --
  ## see sfm.R's identical fix for why.
  call <- match.call()
  model_name <- .match_model_name(model_name, eval(formals()$model_name))
  ## Scale the variance-determinant predictor sits on. ttsfm() has always used
  ## the standard deviation, like sfm(); "var" matches psfm() and the competing
  ## packages, so a delta can be compared across entry points. The default
  ## preserves existing behaviour. See gap C1.
  z_link <- match.arg(z_link)
  .z_sigma <- if (identical(z_link, "sd")) function(eta) exp(eta) else function(eta) sqrt(exp(eta))

  .validate_sfa_call(formula, data, "ttsfm",
    maxit = list(
      maxit.bobyqa = maxit.bobyqa, maxit.psoptim = maxit.psoptim,
      maxit.optim = maxit.optim
    ),
    flags = list(optHessian = optHessian, PSopt = PSopt, inefdec = inefdec)
  )

  .check_model_formula_pipes(formula, model_name)

  DR1 <- data_proc(formula, data, model_name, individual = NULL, inefdec)

  formula <- DR1$formula
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
  if (length(unlist(form_parts)) > 4) { ## might need to incorporate this above
    formula_zp <- DR1$formula_zp
    intercept_zp <- DR1$intercept_zp
    n_zp_vars <- DR1$n_zp_vars
    zp_vars <- DR1$zp_vars
    zp_vars_vec <- DR1$zp_vars_vec
    zp_zp_vec <- DR1$zp_zp_vec
  }


  ## Default starting values for variance equations
  delta <- rep(0.1, length(z_vars))
  delta_p <- rep(0.1, length(zp_vars))
  plm_lm <- lm(formula_x, data_orig)
  beta_0_st <- if (isTRUE(intercept == 0)) {
    NA
  } else {
    plm_lm$coefficients[c(1)]
  }
  beta_hat <- if (isTRUE(intercept == 0)) {
    plm_lm$coefficients[x_vars_vec]
  } else {
    plm_lm$coefficients[x_vars_vec][-1]
  }
  beta_0 <- beta_0_st
  sigma_v <- .2

  ## Starting vector
  if (isTRUE(is.numeric(start_val))) {
    start_v <- start_val
  } else {
    start_v <- if (is.na(beta_0_st)) {
      unname(c(beta_hat, sigma_v, delta, delta_p))
    } else {
      unname(c(beta_0, beta_hat, sigma_v, delta, delta_p))
    }
  }

  ## Output label matrix
  out <- matrix(0, nrow = 3, ncol = length(start_v))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- c(x_vars_vec, "sigv", z_vars, zp_vars)

  DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num = NA, individual = NA, N, model_name, rand.gtre = NULL)

  data <- DR2$data
  Y <- DR2$Y

  data_i_vars <- DR2$data_i_vars
  data_z_vars <- as.matrix(data.frame(subset(data, select = z_vars)))
  data_zp_vars <- as.matrix(data.frame(subset(data, select = zp_vars)))


  if (model_name == "TTNE") {
    fn <- function(p) {
      nr <- n_x_vars ## number of regressors in regression
      nzu <- n_z_vars ## number of determinants for u component
      nzw <- n_zp_vars ## number of determinants for w component

      sigv <- exp(p[nr + 1]) ## Assume homoscedastic two sided component
      sigu <- .z_sigma((data_z_vars %*% p[(nr + 2):(nr + nzu + 1)]))
      sigw <- .z_sigma((data_zp_vars %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)]))

      # if (sigv<= 1e-6){stop("Variance too small")}

      ## Numerical safety (added to fix a real "non-finite value supplied by
      ## optim" error reported from the TTNE example).
      sigu <- pmax(sigu, .SFA_CONSTANTS$MIN_POSITIVE)
      sigw <- pmax(sigw, .SFA_CONSTANTS$MIN_POSITIVE)

      e <- Y - data_i_vars %*% p[1:nr]
      a <- (sigv^2) / (2 * sigw^2) - e / sigw
      b <- e / sigv - sigv / sigw

      alpha <- e / sigu + (sigv^2) / (2 * sigu^2)
      beta <- -e / sigv - sigv / sigu

      ## Clip the exp() arguments before exponentiating.
      alpha[alpha > .SFA_CONSTANTS$EXP_CLIP_UPPER] <- .SFA_CONSTANTS$EXP_CLIP_UPPER
      a[a > .SFA_CONSTANTS$EXP_CLIP_UPPER] <- .SFA_CONSTANTS$EXP_CLIP_UPPER

      denom <- sigu + sigw

      term1 <- exp(alpha)
      term2 <- exp(a)

      ## return will send the summation of the log of the
      ## density of the composed error

      ll <- -log(denom) + log((pnorm(beta) * term1) + (pnorm(b) * term2))

      ## NOTE: fn is passed to minimizers (bobyqa/psoptim/optim all minimize
      ## by default, see opts.R -- none of them flip the sign).
      if (any(is.na(ll))) {
        return(.Machine$double.xmax)
      }
      if (is.null(ll)) {
        return(.Machine$double.xmax)
      }

      ll[ll == -Inf] <- -sqrt(.Machine$double.xmax / length(ll))
      ll[ll == Inf] <- -sqrt(.Machine$double.xmax / length(ll))
      ll[is.nan(ll)] <- -sqrt(.Machine$double.xmax / length(ll))

      return(-sum(ll))
    }

    Start.Time <- start.time()

    prep <- list(n_x_vars, n_z_vars, n_zp_vars)
    names(prep) <- c("n_x_vars", "n_z_vars", "n_zp_vars")

    ## Bounds: .generate_sfa_bounds() returns one lower bound per beta/delta/
    ## delta_p slot.
    lower.BOB0 <- .generate_sfa_bounds(formula, prep)[-c(1:2)]
    lower.BOB <- append(lower.BOB0, -Inf, after = n_x_vars)


    ## ---- Stage 1: bobyqa
    Opt.Bobyqa <- opt.bobyqa(
      fn = fn,
      start_v = start_v,
      lower.bobyqa = lower.BOB,
      maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE,
      verbose = verbose
    )

    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    ## ---- Stage 2: psoptim Every OTHER parameter's window here is
    ## [min(start_v of all other slots) - differ, start_v[j] + differ].
    differ <- 10
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Psoptim <- opt.psoptim(
      fn = fn,
      start_v,
      lower.psoptim = lower1,
      upper.psoptim = c(start_v + differ),
      rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim,
      psopt.TF = PSopt,
      verbose = verbose
    )

    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    ## ---- Stage 3: optim
    differ <- 1
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Optim <- opt.optim(
      fn = fn,
      start_v = start_v,
      lower.optim = lower1,
      upper.optim = c(start_v + differ),
      maxit.optim = maxit.optim,
      opt.TF = optHessian,
      method = Method,
      optHessian = optHessian,
      verbose = verbose
    )

    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    ## Preserve current fallback logic
    if (optHessian == FALSE && PSopt == FALSE) {
      opt <- bob1
    }

    if (optHessian == FALSE && PSopt == TRUE) {
      opt <- opt00
    }

    ## Guard against the "silently accept a failed optim() call" bug described
    ## above.
    if (optHessian == TRUE && !is.null(opt$convergence) && opt$convergence != 0) {
      stop(sprintf(
        "ttsfm() %s: final optimizer stage failed (optim() message: \"%s\"). This can happen when a fit approaches a degenerate boundary (e.g. one variance component -> 0); try a different formula/starting values, or refit with a different random seed if using simulated data.",
        model_name, if (!is.null(opt$message)) opt$message else "unknown"
      ), call. = FALSE)
    }


    ## now for st errs
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
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## metrics
    metrics.ne <- function(p, e = NULL, y = Y, xx = data_i_vars, zu = data_z_vars, zw = data_zp_vars, alphahat = NULL) {
      if (is.null(e)) {
        nr <- ncol(xx) ## Calculate number of regressors in regression
        nzu <- ncol(zu) ## Calculate number of determinants for u component
        nzw <- ncol(zw) ## Calculate number of determinants for w component

        ep.hat <- y - xx %*% p[1:nr]

        if (!is.null(alphahat)) {
          ## Was `alpha.hat` (an undefined global -- flagged by R CMD check's
          ## "no visible binding for global variable" NOTE).
          ep.hat <- y - xx %*% p[1:nr] - alphahat
        }


        sig.v <- exp(p[nr + 1]) ## Assume homoscedastic two sided component
        sig.u <- .z_sigma((zu %*% p[(nr + 2):(nr + nzu + 1)]))
        sig.w <- .z_sigma((zw %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)]))
      } else {
        ep.hat <- e

        sig.v <- exp(p[1])
        sig.u <- exp(p[2])
        sig.w <- exp(p[3])

        ## Have to correct for the shift in our residuals to begin with
        ep.hat <- e - sig.u + sig.w
      }
      ## Use 8.23 and 8.26 to construct metrics
      ## Setup necessary parameters needed.
      lambda <- 1 / sig.w + 1 / sig.u
      a1 <- sig.v^2 / (2 * sig.u^2) + ep.hat / sig.u
      b1 <- -(ep.hat / sig.v + sig.v / sig.u)
      a2 <- sig.v^2 / (2 * sig.w^2) - ep.hat / sig.w
      b2 <- ep.hat / sig.v - sig.v / sig.w
      chi1 <- pnorm(b2) + exp(a1 - a2) * pnorm(b1)
      chi2 <- exp(a2 - a1) * chi1

      Eew.cond <- (lambda / (chi2 * (lambda - 1))) * (pnorm(b1) +
        exp(0.5 * ((b2 + sig.v)^2 - b1^2)) * pnorm(b2 + sig.v))

      Eemw.cond <- (lambda / (chi2 * (1 + lambda))) * (pnorm(b1) +
        exp(a2 - a1 - b2 * sig.v + 0.5 * sig.v^2) * pnorm(b2 - sig.v))

      Eeu.cond <- (lambda / (chi1 * (lambda - 1))) * (pnorm(b2) +
        exp(0.5 * ((b1 + sig.v)^2 - b2^2)) * pnorm(b1 + sig.v))

      Eemu.cond <- (lambda / (chi1 * (1 + lambda))) * (pnorm(b2) +
        exp(a1 - a2 - b1 * sig.v + 0.5 * sig.v) * pnorm(b1 - sig.v))

      Eewmu.cond <- (exp((1 + sig.u) * (a1 + sig.v^2 / 2 / sig.u)) * pnorm(b1 - sig.v) +
        exp((1 - sig.w) * (a2 - sig.v^2 / 2 / sig.w)) *
          pnorm(b2 + sig.v)) / (exp(a1) * pnorm(b1) +
        exp(a2) * pnorm(b2))

      Eeumw.cond <- (exp((1 - sig.u) * (a1 - sig.v^2 / 2 / sig.u)) * pnorm(b1 + sig.v) +
        exp((1 + sig.w) * (a2 + sig.v^2 / 2 / sig.w)) *
          pnorm(b2 - sig.v)) / (exp(a1) * pnorm(b1) +
        exp(a2) * pnorm(b2))

      ## Now calculate the M1 and M2 metrics (these are information deficiency
      ## relative to the actual price) and M5 and M6 metrics.
      M1.ne <- 1 - Eemw.cond
      M2.ne <- Eeu.cond - 1

      M5.ne <- Eew.cond - 1
      M6.ne <- 1 - Eemu.cond

      M7.ne <- Eewmu.cond - 1
      M10.ne <- 1 - Eeumw.cond

      return(list(
        M1 = M1.ne, M2 = M2.ne, M5 = M5.ne, M6 = M6.ne, M7 = M7.ne, M10 = M10.ne, Eew.cond = Eew.cond, Eemw.cond = Eemw.cond,
        Eeu.cond = Eeu.cond, Eemu.cond = Eemu.cond, Eewmu.cond = Eewmu.cond, Eeumw.cond = Eeumw.cond
      ))
    }

    metric.ne.res <- metrics.ne(p = opt$par)

    results <- list(t(out), c(opt), End.Time, start_v, model_name, formula, out["par", ], out["st_err", ], out["t-val", ], metric.ne.res, call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "total_time", "start_v", "model_name", "formula", "coefficients", "std.errors", "t.values", "metrics", "call")
    return(results)
  } else if (model_name == "TTHN") {
    ## Normal - Half Normal - Half Normal two-tier stochastic frontier.
    fn <- function(p) {
      nr <- n_x_vars ## number of regressors in regression
      nzu <- n_z_vars ## number of determinants for u component
      nzw <- n_zp_vars ## number of determinants for w component

      sigv <- exp(p[nr + 1]) ## Assume homoscedastic two sided component
      sigu <- .z_sigma((data_z_vars %*% p[(nr + 2):(nr + nzu + 1)]))
      sigw <- .z_sigma((data_zp_vars %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)]))

      ## Numerical safety, same rationale as the analogous fix in the TTNE
      ## branch above: theta1/theta2/omega1/omega2 below divide by sigv.
      sigv <- pmax(sigv, .SFA_CONSTANTS$MIN_POSITIVE)
      sigu <- pmax(sigu, .SFA_CONSTANTS$MIN_POSITIVE)
      sigw <- pmax(sigw, .SFA_CONSTANTS$MIN_POSITIVE)

      e <- Y - data_i_vars %*% p[1:nr]

      theta1 <- sigw / sigv
      theta2 <- sigu / sigv
      s <- sqrt(sigv^2 + sigw^2 + sigu^2)
      omega1 <- s * sqrt(1 + theta2^2) / theta1
      omega2 <- s * sqrt(1 + theta1^2) / theta2
      lambda1 <- (theta2 / theta1) * sqrt(1 + theta1^2 + theta2^2)
      lambda2 <- (theta1 / theta2) * sqrt(1 + theta1^2 + theta2^2)

      rho1 <- lambda1 / sqrt(1 + lambda1^2)
      rho2 <- -lambda2 / sqrt(1 + lambda2^2)

      x1 <- e / omega1
      x2 <- e / omega2

      ## Bivariate standard normal CDF Phi2(x, 0; rho).
      biv_cdf <- function(xvec, rhovec) {
        mapply(function(xx, rr) {
          mnormt::pmnorm(c(xx, 0), mean = c(0, 0), varcov = matrix(c(1, rr, rr, 1), 2, 2))
        }, xvec, rhovec)
      }

      ## Defensive tryCatch: a pathological parameter draw during optimization
      ## (e.g.
      D <- suppressWarnings(tryCatch(
        biv_cdf(x1, rho1) - biv_cdf(x2, rho2),
        error = function(e) NULL
      ))
      if (is.null(D)) {
        return(.Machine$double.xmax)
      }
      D <- pmax(D, .Machine$double.xmin)

      ll <- log(2 * sqrt(2) / sqrt(pi)) - log(s) - (e^2) / (2 * s^2) + log(D)

      ## Same minimizer-sign convention as the TTNE branch above: return the
      ## NEGATIVE summed log-likelihood (bobyqa/psoptim/optim all minimize fn).
      if (any(is.na(ll))) {
        return(.Machine$double.xmax)
      }
      if (is.null(ll)) {
        return(.Machine$double.xmax)
      }

      ll[ll == -Inf] <- -sqrt(.Machine$double.xmax / length(ll))
      ll[ll == Inf] <- -sqrt(.Machine$double.xmax / length(ll))
      ll[is.nan(ll)] <- -sqrt(.Machine$double.xmax / length(ll))

      return(-sum(ll))
    }

    Start.Time <- start.time()

    prep <- list(n_x_vars, n_z_vars, n_zp_vars)
    names(prep) <- c("n_x_vars", "n_z_vars", "n_zp_vars")

    lower.BOB0 <- .generate_sfa_bounds(formula, prep)[-c(1:2)]
    lower.BOB <- append(lower.BOB0, -Inf, after = n_x_vars)


    ## ---- Stage 1: bobyqa
    Opt.Bobyqa <- opt.bobyqa(
      fn = fn,
      start_v = start_v,
      lower.bobyqa = lower.BOB,
      maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE,
      verbose = verbose
    )

    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    ## ---- Stage 2: psoptim
    differ <- 10
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Psoptim <- opt.psoptim(
      fn = fn,
      start_v,
      lower.psoptim = lower1,
      upper.psoptim = c(start_v + differ),
      rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim,
      psopt.TF = PSopt,
      verbose = verbose
    )

    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    ## ---- Stage 3: optim
    differ <- 1
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Optim <- opt.optim(
      fn = fn,
      start_v = start_v,
      lower.optim = lower1,
      upper.optim = c(start_v + differ),
      maxit.optim = maxit.optim,
      opt.TF = optHessian,
      method = Method,
      optHessian = optHessian,
      verbose = verbose
    )

    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    if (optHessian == FALSE && PSopt == FALSE) {
      opt <- bob1
    }

    if (optHessian == FALSE && PSopt == TRUE) {
      opt <- opt00
    }

    ## See the identical guard in the TTNE branch above for the full
    ## explanation.
    if (optHessian == TRUE && !is.null(opt$convergence) && opt$convergence != 0) {
      stop(sprintf(
        "ttsfm() %s: final optimizer stage failed (optim() message: \"%s\"). This can happen when a fit approaches a degenerate boundary (e.g. one variance component -> 0); try a different formula/starting values, or refit with a different random seed if using simulated data.",
        model_name, if (!is.null(opt$message)) opt$message else "unknown"
      ), call. = FALSE)
    }

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
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## Information-deficiency metrics, generalized to a parameter-vector input
    ## the same way metrics.ne() above is.
    metrics.hn <- function(p, e = NULL, y = Y, xx = data_i_vars, zu = data_z_vars, zw = data_zp_vars) {
      nr <- ncol(xx)
      nzu <- ncol(zu)
      nzw <- ncol(zw)

      if (is.null(e)) {
        ep.hat <- y - xx %*% p[1:nr]
      } else {
        ep.hat <- e
      }

      sig.v <- exp(p[nr + 1])
      sig.u <- .z_sigma((zu %*% p[(nr + 2):(nr + nzu + 1)]))
      sig.w <- .z_sigma((zw %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)]))

      theta1 <- sig.w / sig.v
      theta2 <- sig.u / sig.v
      s <- sqrt(sig.v^2 + sig.u^2 + sig.w^2)
      omega1 <- s * sqrt(1 + theta2^2) / theta1
      omega2 <- s * sqrt(1 + theta1^2) / theta2
      lambda1 <- (theta2 / theta1) * sqrt(1 + theta1^2 + theta2^2)
      lambda2 <- (theta1 / theta2) * sqrt(1 + theta1^2 + theta2^2)

      rho1 <- lambda1 / sqrt(1 + lambda1^2)
      rho2 <- -lambda2 / sqrt(1 + lambda2^2)

      ## General bivariate standard normal CDF Phi2(x, y; rho).
      .biv2 <- function(xvec, yvec, rhovec) {
        n <- max(length(xvec), length(yvec), length(rhovec))
        xvec <- rep_len(xvec, n)
        yvec <- rep_len(yvec, n)
        rhovec <- rep_len(rhovec, n)
        out <- suppressWarnings(tryCatch(
          mapply(function(xx, yy, rr) {
            mnormt::pmnorm(c(xx, yy), mean = c(0, 0), varcov = matrix(c(1, rr, rr, 1), 2, 2))
          }, xvec, yvec, rhovec),
          error = function(e) rep(NA_real_, n)
        ))
        out
      }

      Di <- .biv2(ep.hat / omega1, 0, rho1) - .biv2(ep.hat / omega2, 0, rho2)
      F1i <- 2 * .biv2(ep.hat / omega1, 0, rho1)
      F2i <- 2 * .biv2(ep.hat / omega2, 0, rho2)

      s1 <- sqrt(sig.v^2 + sig.w^2)
      s2 <- sqrt(sig.v^2 + sig.u^2)

      omega.w <- sig.w * s2 / s
      omega.u <- sig.u * s1 / s

      rho.uw <- -sig.w * sig.u / s1 / s2

      Eew.cond <- 2 * (F1i - F2i)^(-1) * exp(0.5 * omega.w^2 + (omega.w / omega1) * ep.hat) *
        (pnorm(-(ep.hat - sig.w^2) / omega2) - .biv2(-(ep.hat - sig.w^2) / omega2, -(omega.w + (ep.hat / omega1)), rho.uw))

      Eemw.cond <- 2 * (F1i - F2i)^(-1) * exp(0.5 * omega.w^2 - (omega.w / omega1) * ep.hat) *
        (pnorm(-(ep.hat + sig.w^2) / omega2) - .biv2(-(ep.hat + sig.w^2) / omega2, (omega.w - (ep.hat / omega1)), rho.uw))

      Eeu.cond <- 2 * (F1i - F2i)^(-1) * exp(0.5 * omega.u^2 - (omega.u / omega2) * ep.hat) *
        (pnorm((ep.hat + sig.u^2) / omega1) - .biv2((ep.hat + sig.u^2) / omega1, ((ep.hat / omega2) - omega.u), rho.uw))

      Eemu.cond <- 2 * (F1i - F2i)^(-1) * exp(0.5 * omega.u^2 + (omega.u / omega2) * ep.hat) *
        (pnorm((ep.hat - sig.u^2) / omega1) - .biv2((ep.hat - sig.u^2) / omega1, ((ep.hat / omega2) + omega.u), rho.uw))

      Eewmu.cond <- exp(((sig.w^2 + sig.u^2) / s^2) * (ep.hat + 0.5 * sig.v^2)) *
        (.biv2((ep.hat + sig.v^2) / omega1, 0, rho1) - .biv2((ep.hat + sig.v^2) / omega2, 0, rho2)) / Di

      Eeumw.cond <- exp(((sig.w^2 + sig.u^2) / s^2) * (0.5 * sig.v^2 - ep.hat)) *
        (.biv2((ep.hat - sig.v^2) / omega1, 0, rho1) - .biv2((ep.hat - sig.v^2) / omega2, 0, rho2)) / Di

      M1 <- 1 - Eemw.cond
      M2 <- Eeu.cond - 1
      M5 <- Eew.cond - 1
      M6 <- 1 - Eemu.cond
      M7 <- Eewmu.cond - 1
      M10 <- 1 - Eeumw.cond

      list(
        M1 = M1, M2 = M2, M5 = M5, M6 = M6, M7 = M7, M10 = M10, Eew.cond = Eew.cond, Eemw.cond = Eemw.cond,
        Eeu.cond = Eeu.cond, Eemu.cond = Eemu.cond, Eewmu.cond = Eewmu.cond, Eeumw.cond = Eeumw.cond
      )
    }

    metric.hn.res <- tryCatch(metrics.hn(p = opt$par), error = function(e) NULL)

    results <- list(t(out), c(opt), End.Time, start_v, model_name, formula, out["par", ], out["st_err", ], out["t-val", ], metric.hn.res, call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "total_time", "start_v", "model_name", "formula", "coefficients", "std.errors", "t.values", "metrics", "call")
    return(results)
  } else if (model_name == "TTNLS") {
    ## Two-tier stochastic frontier via nonlinear least squares.
    fn <- function(p) {
      nr <- n_x_vars
      nzu <- n_z_vars
      nzw <- n_zp_vars

      sigu <- .z_sigma(data_z_vars %*% p[(nr + 2):(nr + nzu + 1)])
      sigw <- .z_sigma(data_zp_vars %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)])

      e <- Y - data_i_vars %*% p[1:nr] + sigu - sigw
      ss <- e^2

      if (any(is.na(ss))) {
        return(.Machine$double.xmax)
      }
      if (is.null(ss)) {
        return(.Machine$double.xmax)
      }

      ss[is.infinite(ss)] <- sqrt(.Machine$double.xmax / length(ss))
      ss[is.nan(ss)] <- sqrt(.Machine$double.xmax / length(ss))

      return(sum(ss))
    }

    Start.Time <- start.time()

    prep <- list(n_x_vars, n_z_vars, n_zp_vars)
    names(prep) <- c("n_x_vars", "n_z_vars", "n_zp_vars")

    lower.BOB0 <- .generate_sfa_bounds(formula, prep)[-c(1:2)]
    lower.BOB <- append(lower.BOB0, -Inf, after = n_x_vars)

    ## ---- Stage 1: bobyqa
    Opt.Bobyqa <- opt.bobyqa(
      fn = fn,
      start_v = start_v,
      lower.bobyqa = lower.BOB,
      maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE,
      verbose = verbose
    )

    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    ## ---- Stage 2: psoptim
    differ <- 10
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Psoptim <- opt.psoptim(
      fn = fn,
      start_v,
      lower.psoptim = lower1,
      upper.psoptim = c(start_v + differ),
      rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim,
      psopt.TF = PSopt,
      verbose = verbose
    )

    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    ## ---- Stage 3: optim
    differ <- 1
    lower1_0 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(n_x_vars + 1)]) - differ)[-c(1:2)]
    lower1 <- append(lower1_0, start_v[n_x_vars + 1] - differ, after = n_x_vars)

    Opt.Optim <- opt.optim(
      fn = fn,
      start_v = start_v,
      lower.optim = lower1,
      upper.optim = c(start_v + differ),
      maxit.optim = maxit.optim,
      opt.TF = optHessian,
      method = Method,
      optHessian = optHessian,
      verbose = verbose
    )

    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    if (optHessian == FALSE && PSopt == FALSE) {
      opt <- bob1
    }

    if (optHessian == FALSE && PSopt == TRUE) {
      opt <- opt00
    }

    ## See the identical guard in the TTNE branch above for the full
    ## explanation.
    if (optHessian == TRUE && !is.null(opt$convergence) && opt$convergence != 0) {
      stop(sprintf(
        "ttsfm() %s: final optimizer stage failed (optim() message: \"%s\"). This can happen when a fit approaches a degenerate boundary (e.g. one variance component -> 0); try a different formula/starting values, or refit with a different random seed if using simulated data.",
        model_name, if (!is.null(opt$message)) opt$message else "unknown"
      ), call. = FALSE)
    }

    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- bob1
      st_err <- rep(NA, length(opt$par))
    }

    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
      st_err <- rep(NA, length(opt$par))
    }

    if (optHessian == TRUE) {
      st_err <- rep(NA_real_, length(opt$par))
      if (!isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
        ## Drop the inert sigv row/column (position n_x_vars+1) before
        ## inverting -- see note above.
        drop_idx <- n_x_vars + 1
        H_sub <- opt$hessian[-drop_idx, -drop_idx, drop = FALSE]
        se_sub <- tryCatch(suppressWarnings(sqrt(diag(solve(H_sub)))),
          error = function(e) rep(NA_real_, nrow(H_sub))
        )
        st_err[-drop_idx] <- se_sub
      }
    }
    ## The scale parameters are NOT IDENTIFIED by this objective, so they
    ## are reported as NA rather than as numbers.
    nls_unident <- (n_x_vars + 1):length(opt$par)
    if (any(is.finite(opt$par[nls_unident]))) {
      warning("ttsfm(model_name = \"TTNLS\"): nonlinear least squares identifies the ",
        "frontier slopes and the composite (intercept + sigma_w - sigma_u) only. ",
        "sigma_v, sigma_u and sigma_w are not separately identified by the ",
        "sum-of-squares objective and are returned as NA; the reported intercept ",
        "absorbs sigma_w - sigma_u. Use model_name = \"TTNE\" or \"TTHN\" if the ",
        "individual scales are needed.",
        call. = FALSE
      )
    }
    opt$par[nls_unident] <- NA_real_
    st_err[nls_unident] <- NA_real_

    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## Information-deficiency metrics: unlike TTNE/TTHN these are direct
    ## functions of sigma_u/sigma_w only, not conditional on epsilon.
    metrics.nls <- function(p, zu = data_z_vars, zw = data_zp_vars) {
      nr <- ncol(data_i_vars)
      nzu <- ncol(zu)
      nzw <- ncol(zw)
      sig.u <- .z_sigma(zu %*% p[(nr + 2):(nr + nzu + 1)])
      sig.w <- .z_sigma(zw %*% p[(nr + nzu + 2):(nr + nzu + nzw + 1)])

      list(
        M1  = 1 - exp(-sig.w),
        M2  = exp(sig.u) - 1,
        M5  = exp(sig.w) - 1,
        M6  = 1 - exp(-sig.u),
        M7  = exp(sig.w - sig.u) - 1,
        M10 = 1 - exp(sig.u - sig.w)
      )
    }

    metric.nls.res <- tryCatch(metrics.nls(p = opt$par), error = function(e) NULL)

    results <- list(t(out), c(opt), End.Time, start_v, model_name, formula, out["par", ], out["st_err", ], out["t-val", ], metric.nls.res, call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "total_time", "start_v", "model_name", "formula", "coefficients", "std.errors", "t.values", "metrics", "call")
    return(results)
  } else {
    stop(paste0(
      "model_name '", model_name, "' is a recognized choice for ttsfm() but has no implementation branch. ",
      "Valid, implemented options are: 'TTNE', 'TTHN', 'TTNLS'."
    ), call. = FALSE)
  }
}
