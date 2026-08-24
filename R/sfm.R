sfm <- function(formula,
                model_name = c(
                  "NHN", "NHN_Z", "NE", "NE_Z", "NR", "THT", "NTN", "NG", "NNAK",
                  "NU", "NGE", "NLN", "NW", "tHN", "TSL"
                ),
                data,
                maxit.bobyqa = 10000,
                maxit.nlminb = 500,
                maxit.psoptim = 1000,
                maxit.optim = 1000,
                REPORT = 1,
                trace = 2,
                pgtol = 0,
                start_val = FALSE,
                PSopt = FALSE,
                use.nlminb = "auto",
                use.bobyqa = "auto",
                optHessian = TRUE,
                inefdec = TRUE,
                upper = NA,
                Method = "L-BFGS-B",
                robust = c("mle", "mlqe", "psi", "mdpd"),
                c_mlqe = 0.20,
                eta = 0.01,
                alpha = 0.2,
                verbose = FALSE,
                Nsim = "auto",
                sim_type = c("halton", "sobol", "torus", "uniform"),
                antithetics = FALSE,
                sim_burn = NULL,
                sim_scrambling = 0L,
                sim_prime = NULL,
                sim_seed = NULL,
                rand.psoptim = NULL,
                keep_objective = FALSE,
                estimator = c("mle", "cols", "mols"),
                cols_boot = 0,
                rand.cols = NULL) {
  ## call/model_name resolution moved ahead of .check_model_formula_pipes().
  call <- match.call()
  model_name <- .match_model_name(model_name, eval(formals()$model_name))
  robust <- match.arg(robust)
  estimator <- match.arg(estimator)
  ## Defaults reproduce the previous behaviour exactly -- Halton, no
  ## antithetics, 1000 discarded -- so no existing result moves.
  sim_type <- match.arg(sim_type)
  if (is.null(sim_burn)) sim_burn <- .SFA_CONSTANTS$HALTON_DISCARD
  if (!is.numeric(sim_burn) || length(sim_burn) != 1L || sim_burn < 0) {
    stop("`sim_burn` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.logical(antithetics) || length(antithetics) != 1L || is.na(antithetics)) {
    stop("`antithetics` must be TRUE or FALSE.", call. = FALSE)
  }

  ## "mols" and "cols" select the same estimator, and the literature's name
  ## for it is MOLS.
  estimator <- if (estimator == "mols") "cols" else estimator

  if (estimator == "cols" && robust != "mle") {
    stop("`robust` applies to the maximum-likelihood estimator only. ",
      "estimator = \"cols\" is a moment estimator and has no divergence to ",
      "robustify; call it with robust = \"mle\".",
      call. = FALSE
    )
  }

  .validate_sfa_call(formula, data, "sfm",
    maxit = list(
      maxit.bobyqa = maxit.bobyqa, maxit.nlminb = maxit.nlminb,
      maxit.psoptim = maxit.psoptim, maxit.optim = maxit.optim
    ),
    flags = list(
      optHessian = optHessian, PSopt = PSopt,
      inefdec = inefdec, verbose = verbose
    )
  )

  .check_model_formula_pipes(formula, model_name)

  ## Robust divergence estimation (MLqE/Psi/MDPD, see R/robust_divergence.R)
  ## is currently only wired up for model_name == "NHN".
  if (robust != "mle" && model_name != "NHN") {
    stop("robust = '", robust, "' is currently only implemented for model_name = ",
      "\"NHN\" (see R/robust_divergence.R). Use model_name = \"NHN\", or ",
      "robust = \"mle\" for other models.",
      call. = FALSE
    )
  }

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

  ## Corrected ordinary least squares.
  if (estimator == "cols") {
    Start.Time <- start.time()
    Xc <- as.matrix(data_i_vars)
    Yc <- inefdec_n * as.numeric(Y)
    ## data_i_vars carries make.names()-mangled labels ("X.Intercept."), while
    ## x_vars_vec keeps the real ones.
    icol <- match("(Intercept)", x_vars_vec)
    if (is.na(icol)) {
      .const <- which(apply(Xc, 2, function(z) length(unique(z)) == 1L))
      icol <- if (length(.const)) .const[[1]] else NA_integer_
    }
    CF <- .cols_fit(Yc, Xc, model_name, intercept_col = icol)

    if (isTRUE(CF$wrong_skew)) {
      warning("sfm(estimator = \"cols\"): the OLS residuals are skewed the WRONG ",
        "way (third central moment ", signif(CF$moments[["m3"]], 3), " >= 0). ",
        "A production frontier implies negative skew, so the moment ",
        "equations have no admissible solution and sigma_u is reported as 0 ",
        "-- read this as no evidence of inefficiency in these data rather ",
        "than as an estimate of zero. This is the Type I failure of Olson, ",
        "Schmidt and Waldman (1980) and is common in small samples.",
        call. = FALSE
      )
    }

    par_v <- c(CF$sigma_v, CF$sigma_u, CF$extra, CF$beta)
    se_v <- c(NA_real_, NA_real_, if (is.null(CF$extra)) NULL else NA_real_, CF$se_beta)
    nm_v <- c("sigv", "sigu", names(CF$extra), x_vars_vec)

    ## Optional nonparametric bootstrap.
    boot_mat <- NULL
    if (is.numeric(cols_boot) && length(cols_boot) == 1L && cols_boot >= 1) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      if (!is.null(rand.cols)) set.seed(rand.cols)
      nobs <- length(Yc)
      boot_mat <- matrix(NA_real_, nrow = as.integer(cols_boot), ncol = length(par_v))
      for (bb in seq_len(as.integer(cols_boot))) {
        idx <- sample.int(nobs, nobs, replace = TRUE)
        bfit <- tryCatch(.cols_fit(Yc[idx], Xc[idx, , drop = FALSE], model_name,
          intercept_col = icol
        ), error = function(e) NULL)
        if (!is.null(bfit)) {
          boot_mat[bb, ] <- c(bfit$sigma_v, bfit$sigma_u, bfit$extra, bfit$beta)
        }
      }
      se_v <- apply(boot_mat, 2, function(z) stats::sd(z, na.rm = TRUE))
      colnames(boot_mat) <- nm_v
    }

    out <- matrix(NA_real_, nrow = 3, ncol = length(par_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- nm_v
    out[1, ] <- par_v
    out[2, ] <- se_v
    out[3, ] <- par_v / se_v

    ## Efficiency at the COLS parameters, using the same predictors the ML path
    ## uses. Undefined when sigma_u collapses under wrong skew.
    eps_c <- inefdec_n * CF$residuals
    if (CF$sigma_u > 0) {
      s2u <- CF$sigma_u^2
      s2v <- CF$sigma_v^2
      exp_u_hat <- .te_battese_coelli(
        mu_star = -eps_c * s2u / (s2u + s2v),
        sigma_star = CF$sigma_u * CF$sigma_v / sqrt(s2u + s2v)
      )
    } else {
      exp_u_hat <- rep(NA_real_, length(eps_c))
    }

    End.Time <- end.time(Start.Time)
    results <- list(
      t(out), NULL, End.Time, NULL, model_name, formula, exp_u_hat,
      CF$wrong_skew, CF$moments, boot_mat, "cols",
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "model_name", "formula",
      "exp_u_hat", "wrong_skew", "residual_moments", "cols_boot_draws",
      "estimator", "coefficients", "std.errors", "t.values", "call"
    )
    return(results)
  }

  ## Fixed low-discrepancy draws for the simulated-ML models (NLN, NW).
  FiMat <- HDraw <- NULL
  if (model_name %in% c("NLN", "NW")) {
    n_obs <- length(as.numeric(Y))

    ## How many simulation draws: the count must grow with n, or the
    ## simulation bias does not vanish.
    .auto_nsim <- function(n) {
      if (model_name == "NLN") {
        max(200L, as.integer(ceiling(3 * sqrt(n))))
      } else {
        max(400L, as.integer(ceiling(8 * sqrt(n))))
      }
    }
    Nsim <- if (identical(Nsim, "auto")) {
      .auto_nsim(n_obs)
    } else {
      max(as.integer(Nsim), 10L)
    }

    Nsim_floor <- .auto_nsim(n_obs)
    if (Nsim < Nsim_floor) {
      warning("sfm(model_name = \"", model_name, "\"): Nsim = ", Nsim,
        " is below the ", Nsim_floor, " that this sample size calls for. ",
        "Simulated ML needs the draw count to grow with n; too few draws ",
        "bias every parameter except the frontier slopes. See ?sfm.",
        call. = FALSE
      )
    }
    ## Draw construction moved to .sml_draws() (matrix_utils.R), which is
    ## shared with the other entry points.
    FiMat <- .sml_draws(n_units = n_obs, n_draws = Nsim, dim = 1L,
                        sim_type = sim_type, antithetics = antithetics,
                        burn = sim_burn, scrambling = sim_scrambling,
                        prime = sim_prime, seed = sim_seed,
                        clamp = 1e-6)[[1L]]
    HDraw <- if (model_name == "NLN") qnorm(FiMat) else -log1p(-FiMat)

    ## NLN integrates in t (see the likelihood), where the integrand is
    ## smooth.
  }

  if (model_name %in% c("NHN", "NE", "NR", "NG", "NNAK", "THT", "NTN", "NHN_Z", "NE_Z", "NU", "NGE", "NLN", "NW", "tHN", "TSL")) {
    like.fn <- function(x) {
      ## The offset is the number of non-beta parameters that precede the
      ## frontier coefficients.
      if (model_name %in% c("NHN", "NE", "NR", "NU", "NGE")) {
        x_x_vec <- x[3:as.numeric(n_x_vars + 2)]
      }
      if (model_name %in% c("THT", "NTN", "NLN", "NW", "tHN", "NG", "NNAK", "TSL")) {
        x_x_vec <- x[4:as.numeric(n_x_vars + 3)]
      }

      if (model_name %in% c("NE_Z", "NHN_Z")) {
        data_z_vars <- as.matrix(data.frame(subset(data, select = z_vars)))
        x_x_vec <- x[2:as.numeric(n_x_vars + 1)]
        z_z_vec <- x[as.numeric(n_x_vars + 2):as.numeric(length(start_v))]
      }

      eps <- (inefdec_n * (Y - as.matrix(data_i_vars) %*% x_x_vec))

      if (model_name == "NHN_Z") {
        sigma_u_fun <- exp(as.matrix(data_z_vars) %*% z_z_vec)
        sigma_v_fun <- x[1]
        sigma_fun <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
        lamb_fun <- sigma_u_fun / sigma_v_fun
        like <- log(pmax((2 / sigma_fun) *
          dnorm(eps / sigma_fun) *
          pnorm(-eps * lamb_fun / sigma_fun), .Machine$double.xmin))
      }

      if (model_name == "NE_Z") {
        sigma_u_fun <- exp(data_z_vars %*% z_z_vec)
        sigv <- x[1]
        l1 <- log(1 / sigma_u_fun)
        l2 <- pnorm(-(eps / sigv) - (sigv / sigma_u_fun), log.p = TRUE)
        l3 <- (eps / sigma_u_fun) + (sigv^2 / (2 * sigma_u_fun^2))
        like <- l1 + l2 + l3
      }

      if (model_name == "NHN") {
        like <- as.numeric(log(pmax((2 / x[2]) *
          dnorm(eps / x[2]) *
          pnorm(-eps * x[1] / x[2]), .Machine$double.xmin)))
      }

      if (model_name == "NE") {
        l1 <- log(1 / x[2])
        l2 <- pnorm(-(eps / x[1]) - (x[1] / x[2]), log.p = TRUE)
        l3 <- (eps / x[2]) + (x[1]^2 / (2 * x[2]^2))
        like <- l1 + l2 + l3
      }

      if (model_name == "NU") {
        ## Normal-uniform (Li 1996, Nguyen 2010): u ~ U(0, theta).
        sigv <- x[1]
        theta <- x[2]
        if (!is.finite(sigv) || !is.finite(theta) || sigv <= 0 || theta <= 0) {
          return(.Machine$double.xmax)
        }
        cdf_hi <- pnorm((eps + theta) / sigv)
        cdf_lo <- pnorm(eps / sigv)
        like <- log(pmax(cdf_hi - cdf_lo, .Machine$double.xmin)) - log(theta)
      }

      if (model_name == "NGE") {
        ## Normal-generalized exponential: u ~ GE(2, lambda), i.e.
        sigv <- x[1]
        sigu <- x[2]
        if (!is.finite(sigv) || !is.finite(sigu) || sigv <= 0 || sigu <= 0) {
          return(.Machine$double.xmax)
        }
        lam <- 1 / sigu
        lt1 <- lam * eps + (lam^2 * sigv^2) / 2 + pnorm(-eps / sigv - lam * sigv, log.p = TRUE)
        lt2 <- 2 * lam * eps + 2 * (lam^2 * sigv^2) + pnorm(-eps / sigv - 2 * lam * sigv, log.p = TRUE)
        d <- pmin(lt2 - lt1, -.Machine$double.eps) ## T2 < T1 by construction
        like <- log(2 * lam) + lt1 + log(-expm1(d))
      }

      if (model_name %in% c("NLN", "NW")) {
        ## Simulated ML.
        sigv <- x[1]
        sigu <- x[2]
        shp <- x[3] ## meanlog (NLN) or shape k (NW)
        ## A large FINITE penalty, not .Machine$double.xmax: optim()'s
        ## finite-difference gradient differences the objective.
        if (!is.finite(sigv) || !is.finite(sigu) || !is.finite(shp) ||
          sigv <= 0 || sigu <= 0 || (model_name == "NW" && shp <= 0)) {
          return(1e12)
        }
        if (model_name == "NLN") {
          ## CHANGE OF VARIABLE, u = sigma_v*t - e.
          e_v <- as.numeric(eps)
          p_hi <- pnorm(e_v / sigv, lower.tail = FALSE)
          Tm <- qnorm(p_hi * (1 - FiMat), lower.tail = FALSE)
          u_draw <- sigv * Tm - e_v ## > 0 by construction
          if (any(!is.finite(u_draw))) {
            return(1e12)
          }
          dens <- p_hi * rowMeans(dlnorm(u_draw, meanlog = shp, sdlog = sigu))
        } else {
          u_draw <- sigu * HDraw^(1 / shp)
          if (any(!is.finite(u_draw))) {
            return(1e12)
          }
          ## eps (length n) recycles down the columns of the n x Nsim matrix.
          dens <- rowMeans(dnorm((as.numeric(eps) + u_draw) / sigv) / sigv)
        }
        like <- log(pmax(dens, .Machine$double.xmin))
      }

      if (model_name == "NR") {
        sigv <- x[1]
        sigu <- x[2]
        sigma <- sqrt(2 * sigv^2 + sigu^2)
        z <- (eps * sigu / sigv) / sigma
        like <- (log(pmax(sigv, .Machine$double.xmin)) - 2 * log(pmax(sigma, .Machine$double.xmin))
          - 1 / 2 * (eps / sigv)^2 + 1 / 2 * z^2 + log(pmax(sqrt(2 / pi) * exp(-1 / 2 * z**2)
            - z * (1 - erf(z / sqrt(2))), .Machine$double.xmin)))
      }

      if (model_name == "THT") {
        sig_u <- x[1]
        sig_v <- x[2]
        a <- x[3]
        lamb <- sig_u / sig_v
        sig <- sqrt(sig_v^2 + sig_u^2)
        ## Tancredi (2002) eq (4): the composed error is skew-t with SCALE
        ## omega = sqrt(sig_v^2 + sig_u^2).
        z_s <- eps / sig
        like <- log(2) - log(pmax(sig, .Machine$double.xmin)) +
          dt(z_s, df = a, log = TRUE) +
          pt(-z_s * lamb * sqrt((a + 1) / (z_s^2 + a)), df = a + 1, log.p = TRUE)
      }

      ## tHN -- Student-t noise with a HALF-NORMAL inefficiency term.
      if (model_name == "tHN") {
        sig_v <- x[1]
        sig_u <- x[2]
        nu <- x[3]
        like <- if (!is.finite(sig_v) || !is.finite(sig_u) || !is.finite(nu) ||
          sig_v <= 0 || sig_u <= 0 || nu <= 2) {
          rep(-.Machine$double.xmax / length(eps), length(eps))
        } else {
          .log_d_thn(eps, sig_v, sig_u, nu)
        }
      }

      if (model_name == "NTN") {
        lam <- x[1]
        sig <- x[2]
        mu <- x[3]

        l1 <- -log(sig^2) / 2
        l2 <- -log(2 * pi) / 2
        l3 <- -(1 / (2 * sig^2)) * (-eps - mu)^2
        l4 <- pnorm(((mu / lam) - eps * lam) / sig, log.p = TRUE)
        l5 <- -pnorm((mu / sig) * sqrt(1 + lam^(-2)), log.p = TRUE)
        like <- l1 + l2 + l3 + l4 + l5
      }

      if (model_name == "TSL") {
        ## Normal / truncated skew-Laplace (Wang 2012). u has density
        sig_v <- x[1]
        sig_u <- x[2]
        lam <- x[3]
        A <- sig_v^2 / (2 * sig_u^2) + eps / sig_u
        B <- (1 + lam)^2 * sig_v^2 / (2 * sig_u^2) + eps * (1 + lam) / sig_u
        aa <- -sig_v / sig_u - eps / sig_v
        bb <- -sig_v * (1 + lam) / sig_u - eps / sig_v
        l1 <- log(2) + A + pnorm(aa, log.p = TRUE)
        l2 <- B + pnorm(bb, log.p = TRUE)
        like <- log1p(lam) - log(2 * lam + 1) - log(sig_u) +
          l1 + log(-expm1(pmin(l2 - l1, -.Machine$double.eps)))
      }

      if (model_name %in% c("NG", "NNAK")) {
        ## Stable log parabolic cylinder function.
        lnDv <- function(nu, z) .log_pcf(nu, z)
        sig_v <- x[1]
        sig_u <- x[2]
        mu <- x[3]
        if (model_name == "NG") {
          like <- ((mu - 1) * log(sig_v) - 1 / 2 * log(2) - 1 / 2 * log(pi) - mu * log(sig_u)
            - 1 / 2 * (eps / sig_v)^2 + 1 / 4 * (eps / sig_v + sig_v / sig_u)^2
            + lnDv(-mu, eps / sig_v + sig_v / sig_u))
        }
        if (model_name == "NNAK") {
          sigma <- sqrt(2 * mu * sig_v^2 + sig_u^2)
          like <- (lgamma(2 * mu) - lgamma(mu) + 1 / 2 * log(2) - 1 / 2 * log(pi) + mu * log(mu)
            + (2 * mu - 1) * log(sig_v) - 2 * mu * log(sigma) - 1 / 2 * (eps / sig_v)^2
            + 1 / 4 * ((eps * sig_u / sig_v) / sigma)^2 + lnDv(-2 * mu, (eps * sig_u / sig_v) / sigma))
        }
      }

      like[like == -Inf] <- -sqrt(.Machine$double.xmax / length(like))
      like[like == Inf] <- -sqrt(.Machine$double.xmax / length(like))
      like[is.nan(like)] <- -sqrt(.Machine$double.xmax / length(like))

      ## Robust divergence estimation (MLqE/Psi/MDPD): swap the standard MLE
      ## objective for the robust one, at these SAME current parameter values.
      if (model_name == "NHN" && robust != "mle") {
        lambda_x <- x[1]
        sigma_x <- x[2]
        sigma_u_x <- (lambda_x * sigma_x) / sqrt(1 + lambda_x^2)
        sigma_v_x <- sigma_u_x / lambda_x
        c_x <- switch(robust,
          mlqe = c_mlqe,
          psi = eta,
          mdpd = alpha
        )
        return(.robust_objective(
          method = robust, loglik = like, c = c_x,
          power_integral_fn = function(p) .nhn_power_integral(sigma_v_x, sigma_u_x, p)
        ))
      }

      return(-sum(like[is.finite(like)]))
    }

    Start.Time <- start.time()

    ## Stage 1: nlminb.
    .nlminb_safe <- c("NHN", "NE", "NTN", "NU")
    if (identical(use.nlminb, "auto")) use.nlminb <- model_name %in% .nlminb_safe
    if (identical(use.bobyqa, "auto")) use.bobyqa <- !(model_name %in% .nlminb_safe)

    grad_fn <- if (model_name == "NHN" && robust == "mle") {
      function(p) .grad_nhn(p, Y, data_i_vars, inefdec_n)
    } else {
      NULL
    }

    ## tHN: multi-start is mandatory, not a refinement.
    thn_starts <- NULL
    if (model_name == "tHN") {
      .thn_try <- function(sv) {
        sv <- pmax(sv, lower_bob + 1e-8)
        o <- tryCatch(suppressWarnings(
          stats::optim(sv, like.fn,
            method = "L-BFGS-B", lower = lower_bob,
            control = list(maxit = 300)
          )
        ), error = function(e) NULL)
        if (is.null(o) || !is.finite(o$value)) NULL else o
      }
      .cand <- list(
        start_v, # OLS/moment start
        replace(start_v, 2:3, c(start_v[2] * 3.0, 4)), # inefficiency-heavy, heavy tail
        replace(start_v, 2:3, c(start_v[2] * 0.25, 25)), # near-normal, light inefficiency
        replace(start_v, 2:3, c(start_v[2] * 8.0, 2.5))
      ) # extreme one-sided
      .fits <- Filter(Negate(is.null), lapply(.cand, .thn_try))
      if (length(.fits)) {
        .obj <- vapply(.fits, function(o) o$value, numeric(1))
        .best <- .fits[[which.min(.obj)]]
        ## distinct optima, to within a tolerance on the objective
        .ndist <- length(unique(round(.obj - min(.obj), 4)))
        thn_starts <- list(
          n_tried = length(.cand), n_converged = length(.fits),
          objectives = -.obj, n_distinct = .ndist,
          best = -min(.obj)
        )
        if (.ndist > 1) {
          warning("sfm(model_name = \"tHN\"): the ", length(.fits), " starting values that ",
            "converged reached ", .ndist, " distinct optima (log-likelihoods ",
            paste(signif(sort(-.obj, decreasing = TRUE), 8), collapse = ", "),
            "). The tHN likelihood is known to be bimodal in nu. The best is ",
            "returned; inspect $thn_starts and consider profiling over a grid ",
            "of fixed nu rather than reporting a single selected value.",
            call. = FALSE
          )
        }
        start_v <- .best$par
      }
    }

    ## Normal-gamma multi-start.
    ng_starts <- NULL
    if (model_name == "NG" && isFALSE(is.numeric(start_val))) {
      .cand <- .ng_start_candidates(epsilon_hat, beta_0_st, beta_hat)
      .cand <- c(list(start_v), .cand)
      .cand <- lapply(.cand, function(z) pmax(z, lower_bob + 1e-10))
      .obj <- vapply(.cand, function(z) {
        tryCatch(
          {
            v <- like.fn(z)
            if (is.finite(v)) v else Inf
          },
          error = function(e) Inf
        )
      }, numeric(1))
      if (any(is.finite(.obj))) {
        ## Polish the most promising few before choosing: a candidate can start in
        ## the right basin yet score worse than one that starts nearer a corner.
        .ord <- order(.obj)[seq_len(min(3L, sum(is.finite(.obj))))]
        .fits <- Filter(Negate(is.null), lapply(.ord, function(i) {
          tryCatch(
            suppressWarnings(
              stats::optim(.cand[[i]], like.fn,
                method = "L-BFGS-B",
                lower = lower_bob, control = list(maxit = 200)
              )
            ),
            error = function(e) NULL
          )
        }))
        .fits <- Filter(function(o) is.finite(o$value), .fits)
        if (length(.fits)) {
          .fo <- vapply(.fits, function(o) o$value, numeric(1))
          .best <- .fits[[which.min(.fo)]]
          ng_starts <- list(
            n_tried = length(.cand), n_polished = length(.fits),
            loglik_at_start = -.obj, best = -min(.fo)
          )
          start_v <- .best$par
        }
      }
    }

    Opt.Nlminb <- opt.nlminb(
      fn = like.fn, start_v = start_v, lower.nlminb = lower_bob,
      gr = grad_fn, maxit.nlminb = maxit.nlminb,
      nlminb.TF = use.nlminb, verbose = verbose
    )
    start_v <- Opt.Nlminb$start_v
    start_feval <- Opt.Nlminb$start_feval

    Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v, lower.bobyqa = lower_bob, maxit.bobyqa = maxit.bobyqa, bob.TF = use.bobyqa, verbose = verbose)
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    Lower.Start <- lower.start(start_v, model_name, differ = 1)

    Opt.Psoptim <- opt.psoptim(
      fn = like.fn, start_v, lower.psoptim = Lower.Start$lower1, rand.psoptim = rand.psoptim,
      upper.psoptim = Lower.Start$upper1, maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    Lower.Start <- lower.start(start_v, model_name, differ = 0.5)
    Opt.Optim <- opt.optim(
      fn = like.fn, start_v = start_v, lower.optim = Lower.Start$lower1,
      upper.optim = Lower.Start$upper1_open, maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
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

    ## A Hessian that cannot be inverted.
    if (optHessian == TRUE) {
      st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
        rep(NA, length(opt$par))
      } else {
        vc <- tryCatch(suppressWarnings(solve(opt$hessian)), error = function(e) NULL)
        if (is.null(vc)) rep(NA, length(opt$par)) else suppressWarnings(sqrt(pmax(diag(vc), 0)))
      }
    }

    ## Robust divergence estimation (MLqE/Psi/MDPD): the naive Hessian-inverse
    ## SE above assumes the information-matrix equality.
    if (model_name == "NHN" && robust != "mle") {
      c_x <- switch(robust,
        mlqe = c_mlqe,
        psi = eta,
        mdpd = alpha
      )
      x_x_idx <- 3:as.numeric(n_x_vars + 2)
      per_obs_fn <- function(par) {
        lambda_p <- par[1]
        sigma_p <- par[2]
        beta_p <- par[x_x_idx]
        sigma_u_p <- (lambda_p * sigma_p) / sqrt(1 + lambda_p^2)
        sigma_v_p <- sigma_u_p / lambda_p
        eps_p <- as.vector(inefdec_n * (Y - as.matrix(data_i_vars) %*% beta_p))
        loglik_p <- .dens_nhn(eps_p, sigma_v_p, sigma_u_p, log = TRUE)
        .robust_objective_vec(
          method = robust, loglik = loglik_p, c = c_x,
          power_integral_fn = function(p) .nhn_power_integral(sigma_v_p, sigma_u_p, p)
        )
      }
      st_err <- if (optHessian == TRUE && !is.null(opt$hessian)) {
        .sandwich_se_nhn(opt$par, opt$hessian, per_obs_fn)
      } else {
        rep(NA, length(opt$par))
      }
    }

    ## `T`/`F` are ordinary variables and can be reassigned by user code, so
    ## the literals TRUE/FALSE are used throughout.
    t_val <- tryCatch(opt$par / st_err, error = function(e) rep(NA_real_, ncol(out)))
    out[1, ] <- opt$par
    out[2, ] <- tryCatch(st_err, error = function(e) rep(NA_real_, ncol(out)))
    out[3, ] <- tryCatch(t_val, error = function(e) rep(NA_real_, ncol(out)))

    ## The posterior of u given the composed residual, for the models where it
    ## is a truncated normal.
    u_post <- NULL

    ## JLMS TE Measurements
    if (model_name %in% c("NHN")) {
      beta <- opt$par[-c(1:2)]
      lamb <- opt$par[1]
      sig <- opt$par[2]
      sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
      sig_v <- sig_u / lamb
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      sig_star <- sig_u * sig_v / sig
      inner <- (lamb * eps_hat) / sig
      exp_u_hat <- ((1 - pnorm((sig_u * sig_v / sig) + inner)) / pmax((1 - pnorm(inner)), .Machine$double.xmin)) *
        exp((sig_u^2 / sig^2) * (eps_hat + 0.5 * sig_v^2))
      exp_u_hat <- pmax(exp_u_hat, 0)
      exp_u_hat <- pmin(exp_u_hat, 1)
      ## Median TE
      mu <- 0
      sigu2 <- sig_u^2
      sigv2 <- sig_v^2
      sig2 <- sigv2 + sigu2
      mu.star <- (sigv2 * mu - sigu2 * (eps_hat)) / sig2
      sig.star <- sig_u * sig_v / sqrt(sig2)

      med_u_hat <- exp(-mu.star + sig.star * qnorm(0.5 * pnorm(mu.star / sig.star)))
      med_u_hat <- pmax(med_u_hat, 0)
      med_u_hat <- pmin(med_u_hat, 1)
      u_post <- list(mu_star = mu.star, sigma_star = rep_len(sig.star, length(mu.star)))
    }

    if (model_name == "NE") {
      ## Normal-exponential.
      beta <- opt$par[-c(1:2)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      mu_star <- -eps_hat - (sig_v^2) / sig_u
      exp_u_hat <- .te_battese_coelli(mu_star, sig_v)
      u_hat <- .jlms_u(mu_star, sig_v)
      u_post <- list(mu_star = mu_star, sigma_star = rep_len(sig_v, length(mu_star)))
    }

    if (model_name == "NTN") {
      ## Normal-truncated normal, u ~ N+(mu, sigma_u^2).
      lamb <- opt$par[1]
      sig <- opt$par[2]
      mu <- opt$par[3]
      beta <- opt$par[-c(1:3)]
      sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
      sig_v <- sig_u / lamb
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      mu_star <- (sig_v^2 * mu - sig_u^2 * eps_hat) / sig^2
      sig_star <- sig_u * sig_v / sig
      exp_u_hat <- .te_battese_coelli(mu_star, sig_star)
      u_hat <- .jlms_u(mu_star, sig_star)
      u_post <- list(mu_star = mu_star, sigma_star = rep_len(sig_star, length(mu_star)))
    }

    if (model_name == "NU") {
      ## u | e is normal with mean -e and sd sigma_v, DOUBLY truncated to [0,
      ## theta] (the uniform prior contributes only the support restriction).
      beta <- opt$par[-c(1:2)]
      sig_v <- opt$par[1]
      theta <- opt$par[2]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      mu_star <- -eps_hat
      a_lo <- (0 - mu_star) / sig_v
      a_hi <- (theta - mu_star) / sig_v
      ## Both CDF differences are taken in log space: for a very efficient
      ## unit a_lo and a_hi are both far into the same tail.
      log_den <- .log_pnorm_diff(a_lo, a_hi)
      log_num <- .log_pnorm_diff(a_lo + sig_v, a_hi + sig_v)
      exp_u_hat <- pmin(pmax(exp(-mu_star + 0.5 * sig_v^2 + log_num - log_den), 0), 1)
      u_hat <- pmax(mu_star + sig_v * (dnorm(a_lo) - dnorm(a_hi)) / exp(log_den), 0)
    }

    if (model_name == "NGE") {
      ## The GE(2, lambda) density is a difference of two exponential terms.
      beta <- opt$par[-c(1:2)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      lam <- 1 / sig_u
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      lt1 <- lam * eps_hat + (lam^2 * sig_v^2) / 2 + pnorm(-eps_hat / sig_v - lam * sig_v, log.p = TRUE)
      lt2 <- 2 * lam * eps_hat + 2 * (lam^2 * sig_v^2) + pnorm(-eps_hat / sig_v - 2 * lam * sig_v, log.p = TRUE)
      ## SIGNED mixture, not a convex one: the GE density is a DIFFERENCE of
      ## two exponential pieces.
      d <- pmin(lt2 - lt1, -.Machine$double.eps)
      w1 <- -1 / expm1(d) ## = 1/(1 - exp(d)),  d < 0 so w1 > 0
      w2 <- w1 - 1
      mu1 <- -eps_hat - lam * sig_v^2
      mu2 <- -eps_hat - 2 * lam * sig_v^2
      exp_u_hat <- pmin(pmax(w1 * .te_battese_coelli(mu1, sig_v) -
        w2 * .te_battese_coelli(mu2, sig_v), 0), 1)
      u_hat <- pmax(w1 * .jlms_u(mu1, sig_v) - w2 * .jlms_u(mu2, sig_v), 0)
    }

    if (model_name %in% c("NLN", "NW")) {
      ## Simulated Bayes rule over the same fixed draws used in the likelihood:
      ##   E[g(u)|e] = sum_s g(u_s) phi((e+u_s)/sigma_v) / sum_s phi((e+u_s)/sigma_v)
      beta <- opt$par[-c(1:3)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      shp <- opt$par[3]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      ## NLN uses the same change of variable as its likelihood (u = sigma_v*t
      ## - e).
      if (model_name == "NLN") {
        e_v <- as.numeric(eps_hat)
        p_hi <- pnorm(e_v / sig_v, lower.tail = FALSE)
        Tm <- qnorm(p_hi * (1 - FiMat), lower.tail = FALSE)
        u_draw <- sig_v * Tm - e_v
        K <- dlnorm(u_draw, meanlog = shp, sdlog = sig_u)
      } else {
        u_draw <- sig_u * HDraw^(1 / shp)
        K <- dnorm((as.numeric(eps_hat) + u_draw) / sig_v) / sig_v
      }
      den <- pmax(rowSums(K), .Machine$double.xmin)
      exp_u_hat <- pmin(pmax(rowSums(K * exp(-u_draw)) / den, 0), 1)
      u_hat <- pmax(rowSums(K * u_draw) / den, 0)
    }

    if (model_name == "NR") {
      beta <- opt$par[-c(1:2)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      sigma <- sqrt(2 * sig_v^2 + sig_u^2)
      z <- (eps_hat * sig_u / sig_v) / sigma
      exp_u_hat <- (exp(1 / 2 * (z + sig_v * sig_u / sigma)^2 - 1 / 2 * z^2) *
        (exp(-1 / 2 * (z + sig_v * sig_u / sigma)^2) - sqrt(pi / 2) * (z + sig_v * sig_u / sigma) * (1 - erf(1 / sqrt(2) * (z + sig_v * sig_u / sigma)))) /
        (exp(-z^2 / 2) - sqrt(pi / 2) * z * (1 - erf(z / sqrt(2)))))
      exp_u_hat <- pmax(exp_u_hat, 0)
    }

    if (model_name == "tHN") {
      ## Bayes rule over the same quadrature nodes the likelihood uses, so
      ## numerator and denominator share the node set.
      beta <- opt$par[-c(1:3)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      nu <- opt$par[3]
      eps_hat <- as.numeric(inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta))))
      .eff <- .thn_eff(eps_hat, sig_v, sig_u, nu)
      exp_u_hat <- .eff$exp_u_hat
      u_hat <- .eff$u_hat

      ## sigma_u collapsing onto zero is a REAL property of this model on real
      ## data, not an optimizer failure.
      thn_nu_unidentified <- isTRUE(nu > 50)
      if (thn_nu_unidentified) {
        warning("sfm(model_name = \"tHN\"): nu converged to ", signif(nu, 4),
          ", where the Student-t noise is numerically indistinguishable from ",
          "normal. The likelihood is flat in nu here, so the value is not an ",
          "estimate -- report it as \"nu large / tails not detectably heavy\", ",
          "and compare against model_name = \"NHN\", which is its limit. The ",
          "other parameters are unaffected. See ?sfm.",
          call. = FALSE
        )
      }

      thn_sigma_u_at_bound <- isTRUE(sig_u / sqrt(sig_u^2 + sig_v^2) < 1e-3)
      if (thn_sigma_u_at_bound) {
        warning("sfm(model_name = \"tHN\"): sigma_u collapsed to the boundary (",
          signif(sig_u, 3), "), so the fit reports essentially no inefficiency ",
          "(mean predicted efficiency ", signif(mean(exp_u_hat, na.rm = TRUE), 4),
          "). This is a known property of the Student-t--half-normal model, not ",
          "necessarily a numerical failure: heavy-tailed noise can absorb the ",
          "entire one-sided component. Inspect $thn_sigma_u_at_bound and treat ",
          "the efficiency scores as uninformative. See ?sfm.",
          call. = FALSE
        )
      }
    }

    if (model_name == "THT") {
      ## Tancredi (2002) section 2.2.
      beta <- opt$par[-c(1:3)]
      sig_u <- opt$par[1]
      sig_v <- opt$par[2]
      a <- opt$par[3]
      eps_hat <- as.numeric(inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta))))
      om2 <- sig_v^2 + sig_u^2
      mu_star <- -eps_hat * sig_u^2 / om2
      s_star <- sqrt(pmax((a + eps_hat^2 / om2) * sig_v^2 * sig_u^2 / (om2 * (a + 1)), .Machine$double.xmin))
      ## E[u|e] has an exact closed form -- it is the mean of a t_(a+1)
      ## truncated to [c, Inf) with c = -mu*/s*, i.e.
      nu_t <- a + 1
      c_trunc <- -mu_star / s_star
      u_hat <- pmax(mu_star + s_star * ((nu_t + c_trunc^2) / (nu_t - 1)) *
        dt(c_trunc, df = nu_t) /
        pmax(pt(c_trunc, df = nu_t, lower.tail = FALSE), .Machine$double.xmin), 0)

      ## The exp(-u) moments have no closed form, so integrate them on the
      ## probability scale over shared quadrature nodes.
      gl <- .gauss_legendre_01(256L)
      F0 <- pt(c_trunc, df = nu_t)
      P <- outer(F0, gl$nodes, function(f, p) f + p * (1 - f))
      Zq <- pmax(mu_star + s_star * qt(P, df = nu_t), 0)
      W <- matrix(gl$weights, nrow = length(mu_star), ncol = length(gl$weights), byrow = TRUE)
      exp_u_hat <- pmin(pmax(rowSums(W * exp(-Zq)), 0), 1)
      sd_exp_u_hat <- sqrt(pmax(rowSums(W * exp(-2 * Zq)) - exp_u_hat^2, 0))
    }

    if (model_name == "NHN_Z") {
      NX <- n_x_vars + 1
      NZ1 <- n_x_vars + 2
      NZ2 <- n_x_vars + n_z_vars + 1
      beta <- opt$par[c(2:NX)]
      delta <- opt$par[c(NZ1:NZ2)]
      sig_v <- opt$par[1]
      sig_u <- exp((as.matrix(as.matrix(data.frame(subset(data, select = z_vars))))) %*% delta)
      lamb <- sig_u / sig_v
      sig <- sqrt(sig_u^2 + sig_v^2)
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      sig_star <- (sig_u * sig_v) / sig
      inner <- (lamb * eps_hat) / sig
      exp_u_hat <- ((1 - pnorm((sig_u * sig_v / sig) + inner)) / pmax((1 - pnorm(inner)), .Machine$double.xmin)) * exp((sig_u^2 / sig^2) * (eps_hat + 0.5 * sig_v^2))
      exp_u_hat <- pmax(exp_u_hat, 0)
      exp_u_hat <- pmin(exp_u_hat, 1)
      u_post <- list(
        mu_star = as.numeric(-eps_hat * sig_u^2 / sig^2),
        sigma_star = as.numeric(sig_star)
      )
    }

    if (model_name == "TSL") {
      ## Truncated skew-Laplace.
      beta <- opt$par[-c(1:3)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      lam <- opt$par[3]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      r1 <- 1 / sig_u
      r2 <- (1 + lam) / sig_u
      lt1 <- log(2) + r1 * eps_hat + (r1^2 * sig_v^2) / 2 +
        pnorm(-eps_hat / sig_v - r1 * sig_v, log.p = TRUE)
      lt2 <- r2 * eps_hat + (r2^2 * sig_v^2) / 2 +
        pnorm(-eps_hat / sig_v - r2 * sig_v, log.p = TRUE)
      d <- pmin(lt2 - lt1, -.Machine$double.eps)
      w1 <- -1 / expm1(d)
      w2 <- w1 - 1
      mu1 <- -eps_hat - r1 * sig_v^2
      mu2 <- -eps_hat - r2 * sig_v^2
      exp_u_hat <- pmin(pmax(w1 * .te_battese_coelli(mu1, sig_v) -
        w2 * .te_battese_coelli(mu2, sig_v), 0), 1)
      u_hat <- pmax(w1 * .jlms_u(mu1, sig_v) - w2 * .jlms_u(mu2, sig_v), 0)
    }

    if (model_name %in% c("NG", "NNAK")) {
      ## Same stable helper the likelihood uses; the series form here also
      ## carried a stray trailing comma inside log().
      lnDv <- function(nu, z) .log_pcf(nu, z)
      beta <- opt$par[-c(1:3)]
      sig_v <- opt$par[1]
      sig_u <- opt$par[2]
      mu <- opt$par[3]
      eps_hat <- inefdec_n * (Y - rowSums(t(t(data_i_vars) * beta)))
      if (model_name == "NG") {
        z <- eps_hat / sig_v + sig_v / sig_u
        exp_u_hat <- exp(((z + sig_v) / 2)^2) / exp((z / 2)^2) * exp(lnDv(-mu, z + sig_v)) / exp(lnDv(-mu, z))
      }
      if (model_name == "NNAK") {
        sigma <- sqrt(2 * mu * sig_v^2 + sig_u^2)
        z <- (eps_hat * sig_u / sig_v) / sigma
        exp_u_hat <- exp((z / 2 + sig_v * sig_u / (2 * sigma))^2) / exp((z / 2)^2) * exp(lnDv(-2 * mu, z + sig_v * sig_u / sigma)) / exp(lnDv(-2 * mu, z))
      }
      exp_u_hat <- pmax(exp_u_hat, 0)
      exp_u_hat <- pmin(exp_u_hat, 1)
    }

    if (model_name %in% c("NHN")) {
      robust_c_report <- if (robust == "mle") {
        NA_real_
      } else {
        switch(robust,
          mlqe = c_mlqe,
          psi = eta,
          mdpd = alpha
        )
      }
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat, med_u_hat,
        out["par", ], out["st_err", ], out["t-val", ], call, robust, robust_c_report
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat", "med_u_hat",
        "coefficients", "std.errors", "t.values", "call", "robust", "robust_c"
      )
    }

    if (model_name %in% c("NHN_Z", "NR", "NNAK")) {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## NG additionally reports its starting-value search, which is what decides
    ## whether the fit lands at the optimum or in the sigma_u = 0 corner.
    if (model_name == "NG") {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat, ng_starts,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat", "ng_starts",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## NE/NTN additionally return u_hat (the JLMS E[u|e] point predictor)
    ## alongside exp_u_hat.
    if (model_name %in% c("NE", "NTN", "NU", "NGE", "NLN", "NW", "TSL")) {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat, u_hat,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat", "u_hat",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## tHN additionally reports the boundary flag and the multi-start diagnostic.
    if (model_name == "tHN") {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat, u_hat,
        thn_sigma_u_at_bound, thn_nu_unidentified, thn_starts,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat", "u_hat",
        "thn_sigma_u_at_bound", "thn_nu_unidentified", "thn_starts",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## THT additionally reports sd_exp_u_hat, the standard deviation of
    ## exp(-u) given the residual (Tancredi 2002, section 2.2).
    if (model_name == "THT") {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula, exp_u_hat, u_hat, sd_exp_u_hat,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula", "exp_u_hat", "u_hat", "sd_exp_u_hat",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## Fallback for models with no efficiency predictor yet (currently the _Z
    ## variants other than NHN_Z).
    if (!(model_name %in% c("NHN", "NHN_Z", "NR", "NG", "NNAK", "NE", "NTN", "NU", "NGE", "NLN", "NW", "THT", "tHN", "TSL"))) {
      results <- list(
        t(out), c(opt), End.Time, start_v, model_name, formula,
        out["par", ], out["st_err", ], out["t-val", ], call
      )
      class(results) <- "sfareg"
      names(results) <- c(
        "out", "opt", "total_time", "start_v", "model_name", "formula",
        "coefficients", "std.errors", "t.values", "call"
      )
    }

    ## Optionally retain the objective, so sfa_diagnostics() can profile the
    ## likelihood and difference it for a gradient after the fact.
    if (isTRUE(keep_objective) && exists("like.fn", inherits = FALSE)) {
      results$objective <- like.fn
    }

    ## Appended here rather than threaded through each per-model list,
    ## which are built positionally and renamed separately.
    if (!is.null(u_post)) {
      results$u_posterior <- u_post
    }

    ## The variance-determinant block, kept so marginal_effects() can report d
    ## E[u]/d z and d Var[u]/d z without re-deriving the design matrix from.
    if (model_name %in% c("NHN_Z", "NE_Z") && isTRUE(n_z_vars > 0)) {
      Zm <- tryCatch(as.matrix(data.frame(subset(data, select = z_vars))),
                     error = function(e) NULL)
      dl <- tryCatch(opt$par[(n_x_vars + 2):length(opt$par)], error = function(e) NULL)
      if (!is.null(Zm) && !is.null(dl) && ncol(Zm) == length(dl)) {
        names(dl) <- colnames(Zm)
        results$z_spec <- list(
          Z = Zm, delta = dl,
          ## sfm()'s _Z models put the linear predictor on the STANDARD
          ## DEVIATION, sigma_u = exp(z'delta).
          link = "sd",
          family = if (model_name == "NHN_Z") "halfnormal" else "exponential"
        )
      }
    }
    return(results)
  } else {
    stop(paste0(
      "model_name '", model_name, "' is a recognized choice for sfm() but has no implementation branch. ",
      "Valid choices are: \"NHN\", \"NHN_Z\", \"NE\", \"NE_Z\", \"NR\", \"THT\", \"NTN\", \"NG\", \"NNAK\", \"NU\", \"NGE\", \"NLN\", \"NW\", \"tHN\", \"TSL\"."
    ), call. = FALSE)
  }
}
