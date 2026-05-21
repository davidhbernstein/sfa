.gtre_prepare <- function(data,
                          individual,
                          y_var,
                          x_vars_vec,
                          z_vars,
                          zp_vars,
                          n_x_vars,
                          beta_hat,
                          beta_0,
                          beta_0_st,
                          sigma_v,
                          sigma_r,
                          start_val,
                          halton_num,
                          rand.gtre,
                          N,
                          verbose,
                          formula) {
  
  ## Default starting values for variance equations
  delta   <- rep(0.1, length(z_vars))
  delta_p <- rep(0.1, length(zp_vars))
  
  ## Starting vector
  if(isTRUE(is.numeric(start_val))){
    start_v <- start_val
  } else {
    start_v <- if(is.na(beta_0_st)){
      unname(c(sigma_v, sigma_r, beta_hat, delta, delta_p))
    } else {
      unname(c(sigma_v, sigma_r, beta_0, beta_hat, delta, delta_p))
    }
  }
  
  ## Output label matrix
  out           <- matrix(0, nrow = 3, ncol = length(start_v))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- c("sigv", "sigr", colnames(data[, x_vars_vec, drop = FALSE]), z_vars, zp_vars)
  
  ## Number of simulation draws
  R <- if(isTRUE(is.numeric(halton_num))) halton_num else ceiling(sqrt(nrow(data))) + 100
  
  ## Halton draws
  R_H <- randtoolbox::halton(
    R + .SFA_CONSTANTS$HALTON_DISCARD,
    2,
    start = 1,
    normal = FALSE
  )[-c(1:.SFA_CONSTANTS$HALTON_DISCARD), c(1:2)]
  
  ## First column: standard normal via qnorm
  ## Second column: half-normal via inverse error function
  R_H <- cbind(
    qnorm(R_H[,1]),
    sqrt(2) * pracma::erfinv(R_H[,2])
  )
  
  if(!is.null(rand.gtre)){
    set.seed(rand.gtre)
  }
  
  ## Optional decorrelation step retained from current code
  mat <- matrix(0, nrow = R, ncol = 9999)
  for(v in 1:9999){
    mat[,v] <- sample(R_H[,1])
  }
  
  cor_vec <- matrix(0, 9999, 1)
  for(v in 1:9999){
    cor_vec[v] <- abs(cor(mat[,v], R_H[,2]))
  }
  
  R_H <- cbind(mat[, which(cor_vec == min(cor_vec))[1]], R_H[,2])
  rm(mat, cor_vec, v)
  
  ## Build firm-level lists
  indiv <- noquote(as.vector(unique(data[, individual])))
  t_vec <- rep(0, N)
  
  data_i       <- vector("list", N)
  Y            <- vector("list", N)
  data_i_vars  <- vector("list", N)
  data_z_vars  <- vector("list", N)
  data_zp_vars <- vector("list", N)
  R_h1         <- vector("list", N)
  R_h2         <- vector("list", N)
  
  for(ii in seq_len(N)){
    data_i[[ii]]       <- data[which(data[, individual] == indiv[ii]), , drop = FALSE]
    t_vec[ii]          <- nrow(data_i[[ii]])
    R_h1[[ii]]         <- t(matrix(rep(R_H[,1], t_vec[[ii]]), R, t_vec[[ii]]))
    R_h2[[ii]]         <- abs(t(matrix(rep(R_H[,2], t_vec[[ii]]), R, t_vec[[ii]])))
    Y[[ii]]            <- matrix(rep(data_i[[ii]][, y_var], R), t_vec[[ii]], R)
    data_i_vars[[ii]]  <- data.frame(data_i[[ii]][, x_vars_vec, drop = FALSE])
    data_z_vars[[ii]]  <- data.frame(data_i[[ii]][, z_vars, drop = FALSE])
    data_zp_vars[[ii]] <- data.frame(data_i[[ii]][, zp_vars, drop = FALSE])
  }
  
  list(
    data         = data,
    individual   = individual,
    y_var        = y_var,
    x_vars_vec   = x_vars_vec,
    z_vars       = z_vars,
    zp_vars      = zp_vars,
    n_x_vars     = n_x_vars,
    n_z_vars     = length(z_vars),
    n_zp_vars    = length(zp_vars),
    N            = N,
    R            = R,
    start_v      = start_v,
    out_template = out,
    indiv        = indiv,
    t            = t_vec,
    data_i       = data_i,
    Y            = Y,
    data_i_vars  = data_i_vars,
    data_z_vars  = data_z_vars,
    data_zp_vars = data_zp_vars,
    R_h1         = R_h1,
    R_h2         = R_h2,
    formula      = formula 
  )
}


## ------------------------------------------------------------
## Helper: unpack GTRE parameter vector
##
## Assumes ordering:
##   1                        -> sigma_v
##   2                        -> sigma_r
##   3:(2+n_x_vars)           -> beta
##   next n_z_vars            -> delta
##   next n_zp_vars           -> delta_p
## ------------------------------------------------------------
.gtre_unpack_par <- function(x, prep){
  beta_start   <- 3
  beta_end     <- beta_start + prep$n_x_vars - 1
  
  delta_start  <- beta_end + 1
  delta_end    <- delta_start + prep$n_z_vars - 1
  
  deltap_start <- delta_end + 1
  deltap_end   <- deltap_start + prep$n_zp_vars - 1
  
  beta <- x[beta_start:beta_end]
  
  delta <- if(prep$n_z_vars > 0) {
    x[delta_start:delta_end]
  } else {
    numeric(0)
  }
  
  delta_p <- if(prep$n_zp_vars > 0) {
    x[deltap_start:deltap_end]
  } else {
    numeric(0)
  }
  
  list(
    sig_v   = x[1],
    sig_r   = x[2],
    beta    = beta,
    delta   = delta,
    delta_p = delta_p
  )
}


## ------------------------------------------------------------
## GTRE simulated negative log-likelihood
##
## This is the computational core to optimize later.
##
## Important design choice:
##   This function is "pure" with respect to the preparation
##   object prep. That makes it much easier to benchmark and
##   optimize later.
## ------------------------------------------------------------
.gtre_loglik <- function(x, prep, inefdec_n){
  
  par <- .gtre_unpack_par(x, prep)
  
  ## Enforce positivity on sigma_v and sigma_r
  sig_v <- max(par$sig_v, .SFA_CONSTANTS$MIN_POSITIVE)
  sig_r <- max(par$sig_r, .SFA_CONSTANTS$MIN_POSITIVE)
  
  ## Bound the variance linear predictors before exponentiating.
  ## This avoids exp() overflow when the optimizer visits extreme
  ## parameter/Z combinations.  exp(40) is already enormous for
  ## variance-scale calculations, so this is a numerical guard, not
  ## an economically meaningful restriction in ordinary cases.
  eta_bound <- 40
  
  ## One firm-level log-likelihood contribution
  ll_i <- function(ii){
    
    Ti <- prep$t[ii]
    
    ## Persistent inefficiency scale for firm i
    ##
    ## IMPORTANT:
    ##   This is computed inside the firm-specific contribution,
    ##   not outside, so it is truly unit-specific.
    sigma_h_fun <- if(prep$n_zp_vars > 0){
      eta_h_i <- as.numeric(as.matrix(prep$data_zp_vars[[ii]]) %*% par$delta_p)
      if(any(!is.finite(eta_h_i))){
        return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
      }
      eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
      mean(sqrt(exp(eta_h_i)))
    } else {
      1
    }
    
    ## Construct epsilon_it draw-by-draw
    eps_ii <- prep$Y[[ii]] - sig_r * prep$R_h1[[ii]] + sigma_h_fun * prep$R_h2[[ii]] * inefdec_n
    
    ## Remove frontier mean
    for(qq in seq_len(prep$n_x_vars)){
      eps_ii <- eps_ii - par$beta[qq] *
        matrix(rep(prep$data_i_vars[[ii]][, qq], prep$R), Ti, prep$R)
    }
    
    eps_ii <- inefdec_n * eps_ii
    
    ## Transient inefficiency scale
    sigma_u_fun <- if(prep$n_z_vars > 0){
      eta_u_i <- as.numeric(as.matrix(prep$data_z_vars[[ii]]) %*% par$delta)
      if(any(!is.finite(eta_u_i))){
        return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
      }
      eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
      sqrt(exp(eta_u_i))
    } else {
      rep(1, Ti)
    }
    
    sigma_fun <- sqrt(sig_v^2 + sigma_u_fun^2)
    lamb_fun  <- sigma_u_fun / sig_v
    
    ## Expand firm-level sigma objects across draws
    sigma_mat <- matrix(rep(sigma_fun, prep$R), Ti, prep$R)
    lamb_mat  <- matrix(rep(lamb_fun,  prep$R), Ti, prep$R)
    
    ## Simulated likelihood contribution:
    ##   average over R draws of the product across t
    ##
    ## NOTE:
    ##   This is still written close to your current code to
    ##   preserve behavior before optimization/refactoring.
    sim_terms <- (2 / sigma_mat) *
      dnorm(eps_ii / sigma_mat) *
      pmax(
        pnorm(-eps_ii * lamb_mat / sigma_mat),
        eps_ii * 0 + .SFA_CONSTANTS$MIN_POSITIVE
      )
    
    prod_vec_n <- log(mean(colProds(sim_terms)))
    
    -prod_vec_n
  }
  
  ll_vec <- unlist(lapply(seq_len(prep$N), ll_i))
  
  ll_vec[which(ll_vec == Inf)]  <- (.SFA_CONSTANTS$MAX_VALUE)^0.1
  ll_vec[which(ll_vec == -Inf)] <- -(.SFA_CONSTANTS$MAX_VALUE)^0.1
  
  sum(ll_vec[is.finite(ll_vec)])
}


## ------------------------------------------------------------
## Staged optimization for GTRE
##
## Preserves your current bobyqa -> psoptim -> optim sequence.
## ------------------------------------------------------------
.gtre_optimize <- function(prep,
                           inefdec_n,
                           maxit.bobyqa,
                           rand.psoptim,
                           maxit.psoptim,
                           PSopt,
                           maxit.optim,
                           optHessian,
                           Method,
                           verbose){
  
  fn <- function(x) .gtre_loglik(x = x, prep = prep, inefdec_n = inefdec_n)
  
  Start.Time <- start.time()
  
  lower.BOB <- .generate_sfa_bounds(prep$formula, prep)  # default to -Inf 
  
  ## ---- Stage 1: bobyqa
  Opt.Bobyqa <- opt.bobyqa(
    fn = fn,
    start_v = prep$start_v,
    lower.bobyqa = lower.BOB,
    maxit.bobyqa = maxit.bobyqa,
    bob.TF = TRUE,
    verbose = verbose
  )
  
  start_v     <- Opt.Bobyqa$start_v
  start_feval <- Opt.Bobyqa$start_feval
  bob1        <- Opt.Bobyqa$bob1
  
  ## ---- Stage 2: psoptim
  differ <- 10
  # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
  lower1 <- .generate_sfa_bounds(prep$formula, prep , inf_sub = min(start_v[-c(1:2)])-differ )   
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
  
  start_v     <- Opt.Psoptim$start_v
  start_feval <- Opt.Psoptim$start_feval
  opt00       <- Opt.Psoptim$opt00
  
  ## ---- Stage 3: optim
  differ <- 1
  # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
  lower1 <- .generate_sfa_bounds(prep$formula, prep , inf_sub = min(start_v[-c(1:2)])-differ )   
  
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
  
  start_v     <- Opt.Optim$start_v
  start_feval <- Opt.Optim$start_feval
  opt         <- Opt.Optim$opt
  
  End.Time <- end.time(Start.Time)
  
  ## Preserve current fallback logic
  if(optHessian == FALSE && PSopt == FALSE){
    opt <- bob1
  }
  
  if(optHessian == FALSE && PSopt == TRUE){
    opt <- opt00
  }
  
  list(
    fn          = fn,
    opt         = opt,
    start_v     = start_v,
    start_feval = start_feval,
    total_time  = End.Time
  )
}


## ------------------------------------------------------------
## Post-estimation GTRE technical efficiency recovery
##
## This uses the refactored TE code we developed earlier.
## ------------------------------------------------------------
.gtre_te <- function(opt, prep, inefdec_n){
  
  ## Basic dimensions and indexing
  n       <- sum(prep$t)
  id_obs  <- rep(seq_len(prep$N), prep$t)
  t_cum   <- cumsum(prep$t)
  t_start <- c(1, head(t_cum, -1) + 1)
  
  ## Build variance-design matrices
  Z_mat  <- .make_var_design(prep$data, prep$z_vars,  rows = NULL,    int_name = "int_u")
  Zp_mat <- .make_var_design(prep$data, prep$zp_vars, rows = t_start, int_name = "int_h")
  
  # n_z_eff  <- ncol(Z_mat)
  # n_zp_eff <- ncol(Zp_mat)
  
  ## Parameter indexing
  beta_start   <- 3
  beta_end     <- beta_start + prep$n_x_vars - 1
  
  delta_start  <- beta_end + 1
  delta_end    <- delta_start + prep$n_z_vars - 1
  
  deltap_start <- delta_end + 1
  deltap_end   <- deltap_start + prep$n_zp_vars - 1
  
  beta    <- opt$par[beta_start:beta_end]
  delta   <- opt$par[delta_start:delta_end]
  delta_p <- opt$par[deltap_start:deltap_end]
  
  sig_v <- max(opt$par[1], 1e-8)
  sig_r <- max(opt$par[2], 1e-8)
  
  min_sd    <- 1e-8
  eta_bound <- 40
  
  ## Compute and check the variance linear predictors before exp().
  ## The earlier pmax() guard only protected against near-zero standard
  ## deviations.  It did not protect against overflow from exp(Z %*% delta),
  ## which can produce Inf values and cause the posterior covariance system
  ## to fail later in .safe_inverse().
  eta_u <- as.numeric(Z_mat  %*% delta)
  eta_h <- as.numeric(Zp_mat %*% delta_p)
  
  if(any(!is.finite(eta_u)) || any(!is.finite(eta_h))){
    stop(
      "Non-finite variance linear predictor in .gtre_te(). Check Z*delta and Zp*delta_p.",
      call. = FALSE
    )
  }
  
  ## Bound before exponentiating to prevent numerical overflow.
  eta_u <- pmin(pmax(eta_u, -eta_bound), eta_bound)
  eta_h <- pmin(pmax(eta_h, -eta_bound), eta_bound)
  
  sig_u_all <- pmax(sqrt(exp(eta_u)), min_sd)
  sig_h_all <- pmax(sqrt(exp(eta_h)), min_sd)
  
  if(any(!is.finite(sig_u_all)) || any(!is.finite(sig_h_all))){
    stop(
      "Non-finite sigma_u or sigma_h values in .gtre_te(). Check bounded variance predictors.",
      call. = FALSE
    )
  }
  
  sig_u_split <- split(sig_u_all, id_obs)
  
  ## Residual vectors epsilon_i
  e_i <- Map(
    f = function(Yi, Xi){
      pmin(
        inefdec_n * (Yi[,1] - rowSums(t(t(Xi) * beta))),
        Yi[,1] * 0
      )
    },
    Yi = prep$Y,
    Xi = prep$data_i_vars
  )
  
  ## Firm-level matrices
  A_i <- lapply(prep$t, function(Ti) -cbind(rep(1, Ti), diag(Ti)))
  
  SIG <- lapply(
    prep$t,
    function(Ti) sig_v^2 * diag(Ti) + sig_r^2 * tcrossprod(rep(1, Ti))
  )
  
  VEE <- Map(
    f = function(sig_h_i, sig_u_i){
      Ti <- length(sig_u_i)
      rbind(
        c(sig_h_i^2, rep(0, Ti)),
        cbind(rep(0, Ti), diag(sig_u_i^2, nrow = Ti, ncol = Ti))
      )
    },
    sig_h_i = sig_h_all,
    sig_u_i = sig_u_split
  )
  
  ## Build posterior covariance objects one firm at a time so that
  ## any numerical failure reports the offending firm and panel length.
  ## This makes failures much easier to debug than a generic Map/mapply
  ## traceback.
  post_obj <- lapply(seq_len(prep$N), function(ii){
    if(any(!is.finite(VEE[[ii]])) || any(!is.finite(A_i[[ii]])) || any(!is.finite(SIG[[ii]]))){
      stop(
        sprintf(
          "Non-finite input matrix in .gtre_te() for firm %s (Ti = %s).",
          ii, prep$t[ii]
        ),
        call. = FALSE
      )
    }
    
    tryCatch(
      .safe_linear_combo(
        VEE_i      = VEE[[ii]],
        A_i        = A_i[[ii]],
        SIG_i      = SIG[[ii]],
        base_ridge = 1e-10,
        ridge_mult = 10,
        max_tries  = 8
      ),
      error = function(e){
        stop(
          sprintf(
            "Failed to compute GTRE posterior matrices in .gtre_te() for firm %s (Ti = %s): %s",
            ii, prep$t[ii], conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  })
  
  LAM <- lapply(post_obj, `[[`, "LAM")
  ARR <- lapply(post_obj, `[[`, "ARR")
  
  ridge_report <- data.frame(
    firm        = seq_len(prep$N),
    ridge_VEE   = sapply(post_obj, `[[`, "ridge_VEE"),
    ridge_SIG   = sapply(post_obj, `[[`, "ridge_SIG"),
    ridge_K     = sapply(post_obj, `[[`, "ridge_K"),
    method_VEE  = sapply(post_obj, `[[`, "invVEE_method"),
    method_SIG  = sapply(post_obj, `[[`, "invSIG_method"),
    method_K    = sapply(post_obj, `[[`, "invK_method")
  )
  
  ## Persistent TE
  res_d <- mapply(
    FUN = function(Ti, ARR_i, e_i, LAM_i){
      ptmvnorm(
        lowerx = rep(0, Ti + 1),
        upperx = rep(Inf, Ti + 1),
        mean   = as.numeric(ARR_i %*% e_i),
        sigma  = LAM_i
      )[1]
    },
    Ti       = prep$t,
    ARR_i    = ARR,
    e_i      = e_i,
    LAM_i    = LAM,
    SIMPLIFY = FALSE
  )
  
  res_n <- mapply(
    FUN = function(Ti, ARR_i, e_i, LAM_i){
      shift_vec <- c(-1, rep(0, Ti))
      ptmvnorm(
        lowerx = rep(0, Ti + 1),
        upperx = rep(Inf, Ti + 1),
        mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
        sigma  = LAM_i
      )[1]
    },
    Ti       = prep$t,
    ARR_i    = ARR,
    e_i      = e_i,
    LAM_i    = LAM,
    SIMPLIFY = FALSE
  )
  
  H <- mapply(
    FUN = function(Ti, ARR_i, e_i, LAM_i, rd, rn){
      shift_vec <- c(-1, rep(0, Ti))
      (max(rn, .SFA_CONSTANTS$MIN_POSITIVE) /
          max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
        exp(
          t(shift_vec) %*% ARR_i %*% e_i +
            0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
        )
    },
    Ti    = prep$t,
    ARR_i = ARR,
    e_i   = e_i,
    LAM_i = LAM,
    rd    = res_d,
    rn    = res_n
  )
  
  H <- pmin(H, 1)
  
  ## Transient TE
  U_list <- mapply(
    FUN = function(Ti, ARR_i, e_i, LAM_i, rd){
      
      sapply(seq_len(Ti), function(j){
        
        shift_vec <- rep(0, Ti + 1)
        shift_vec[j + 1] <- -1
        
        rn_t <- ptmvnorm(
          lowerx = rep(0, Ti + 1),
          upperx = rep(Inf, Ti + 1),
          mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
          sigma  = LAM_i
        )[1]
        
        (max(rn_t, .SFA_CONSTANTS$MIN_POSITIVE) /
            max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
          exp(
            t(shift_vec) %*% ARR_i %*% e_i +
              0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
          )
      })
    },
    Ti       = prep$t,
    ARR_i    = ARR,
    e_i      = e_i,
    LAM_i    = LAM,
    rd       = res_d,
    SIMPLIFY = FALSE
  )
  
  U <- unlist(U_list, use.names = FALSE)
  U <- pmin(U, 1)
  
  list(
    U = U,
    H = H,
    ridge_report = ridge_report
  )
}


## ------------------------------------------------------------
## Finalize GTRE result object
##
## Handles:
##   - standard errors
##   - t-values
##   - result labeling
##   - sfareg return object
## ------------------------------------------------------------
.gtre_finalize <- function(opt_obj,
                           prep,
                           data,
                           formula,
                           model_name,
                           call,
                           U,
                           H){
  
  opt <- opt_obj$opt
  out <- prep$out_template
  
  ## Standard errors
  if(isTRUE(any(opt$hessian == 0))){
    st_err <- rep(NA, length(opt$par))
  } else if(is.null(opt$hessian)) {
    st_err <- rep(NA, length(opt$par))
  } else {
    st_err <- suppressWarnings(
      sqrt(pmax(diag(solve(opt$hessian)), 0))
    )
  }
  
  t_val   <- opt$par / st_err
  out[1,] <- opt$par
  out[2,] <- st_err
  out[3,] <- t_val
  
  ## Rename intercept labels if needed
  if(length(colnames(out)[which(colnames(out) == "(Intercept)")]) > 2){
    colnames(out)[which(colnames(out) == "(Intercept)")] <-
      c("(Intercept x)", "(Intercept u)", "(Intercept h)")
  }
  
  results <- list(
    t(out),
    c(opt),
    data,
    opt_obj$total_time,
    opt_obj$start_v,
    model_name,
    formula,
    U,
    H,
    out["par",],
    out["st_err",],
    out["t-val",],
    call
  )
  
  class(results) <- "sfareg"
  names(results) <- c(
    "out",
    "opt",
    "data",
    "total_time",
    "start_v",
    "model_name",
    "formula",
    "U",
    "H",
    "coefficients",
    "std.errors",
    "t.values",
    "call"
  )
  
  results
}


## ------------------------------------------------------------
## Main internal GTRE wrapper
##
## This is the single entry point for the GTRE branch.
##
## Recommended usage:
##   results <- .gtre(...)
##
## Then the higher-level package code can dispatch by model_name.
## ------------------------------------------------------------
.gtre <- function(data,
                  individual,
                  y_var,
                  x_vars_vec,
                  z_vars,
                  zp_vars,
                  n_x_vars,
                  beta_hat,
                  beta_0,
                  beta_0_st,
                  sigma_v,
                  sigma_r,
                  start_val,
                  halton_num,
                  rand.gtre,
                  N,
                  inefdec_n,
                  maxit.bobyqa,
                  rand.psoptim,
                  maxit.psoptim,
                  PSopt,
                  maxit.optim,
                  optHessian,
                  Method,
                  formula,
                  call,
                  verbose = FALSE){
  
  ## 1. Static setup
  prep <- .gtre_prepare(
    data       = data,
    individual = individual,
    y_var      = y_var,
    x_vars_vec = x_vars_vec,
    z_vars     = z_vars,
    zp_vars    = zp_vars,
    n_x_vars   = n_x_vars,
    beta_hat   = beta_hat,
    beta_0     = beta_0,
    beta_0_st  = beta_0_st,
    sigma_v    = sigma_v,
    sigma_r    = sigma_r,
    start_val  = start_val,
    halton_num = halton_num,
    rand.gtre  = rand.gtre,
    N          = N,
    verbose    = verbose,
    formula    = formula 
  )
  
  ## 2. Optimization
  opt_obj <- .gtre_optimize(
    prep          = prep,
    inefdec_n     = inefdec_n,
    maxit.bobyqa  = maxit.bobyqa,
    rand.psoptim  = rand.psoptim,
    maxit.psoptim = maxit.psoptim,
    PSopt         = PSopt,
    maxit.optim   = maxit.optim,
    optHessian    = optHessian,
    Method        = Method,
    verbose       = verbose
  )
  
  ## 3. Post-estimation technical efficiency
  te_obj <- .gtre_te(
    opt       = opt_obj$opt,
    prep      = prep,
    inefdec_n = inefdec_n
  )
  
  ## 4. Finalize result object
  .gtre_finalize(
    opt_obj    = opt_obj,
    prep       = prep,
    data       = data,
    formula    = formula,
    model_name = "GTRE_Z",
    call       = call,
    U          = te_obj$U,
    H          = te_obj$H
  )
}




