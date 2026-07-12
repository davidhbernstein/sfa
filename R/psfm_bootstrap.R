## ============================================================================
## psfm_bootstrap()
## Parametric bootstrap for panel stochastic frontier models ( GTRE /
## GTRE with parametric h) fit with psfm() from the sfm package.
##
## Authors:  David H. Bernstein & Chris Parmeter
## Adapted into a reusable, package-ready function: 11/2025
##
## ----------------------------------------------------------------------------
## WHAT THIS DOES (matches the original bootstrap logic exactly, generalized):
##   1. Pull frontier (x), inefficiency-Z (u), and inefficiency-Z (h) blocks
##      directly from the model's three-part Formula (y ~ x | z_u | z_h).
##      The h-part is OPTIONAL: if the fitted formula has only two RHS parts,
##      h is treated either as a single scalar sigma (TRE-style) or omitted
##      entirely, controlled by `h_type`.
##   2. Re-draw v, u (half-normal w/ z-covariates), r (individual effects),
##      and h (half-normal, either parametric-in-covariates or scalar) exactly
##      as the original model assumes. Parameterization confirmed from
##      sfa:::psfm source: sigma_u_i = sqrt(exp(z %*% delta)),
##      sigma_h_i = sqrt(exp(h_vars %*% delta_p)).
##   3. Build a new response y* = X * beta_hat + v + u + r + h
##   4. Re-estimate psfm() on the simulated data, in parallel, BOOT times.
##   5. Return bootstrap SEs / t-values for all parameters in $out, plus the
##      raw bootstrap draws (parameter matrix + efficiency-score matrices).
##
## NOTE ON ROW LAYOUT OF psfm_object$out (confirmed structure):
##   row 1        : sigv
##   row 2        : sigr
##   next K_x rows: frontier (x) parameters figure, in formula-part-1 order
##   next K_z rows: u-block z parameters, in formula-part-2 order
##   next K_h rows: h-block z parameters, in formula-part-3 order (if present)
## ============================================================================

psfm_bootstrap <- function(psfm_object,
                          numCores,
                          BOOT,
                          individual,
                          h_type        = c("auto", "none", "scalar", "parametric"),
                          maxit.psoptim = 1000,
                          seed_offset   = 0,
                          write_back    = TRUE,
                          pkgs          = c("sfa", "Formula", "pbapply"),
                          inefdec,
                          rand.gtre     = NULL,
                          rand.psoptim  = NULL,
                          maxit.bobyqa  = 1,
                          maxit.optim   = 1 ) {
  
  h_type <- match.arg(h_type)
  
  ## ---- 0. Basic validation -------------------------------------------------
  required_fields <- c("out", "data", "formula", "model_name")
  missing_fields  <- setdiff(required_fields, names(psfm_object))
  if (length(missing_fields) > 0) {
    stop("psfm_object is missing required component(s): ",
         paste(missing_fields, collapse = ", "))
  }
  
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required to parse the multi-part model formula.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required to run the bootstrap in parallel.")
  }
  
  data   <- psfm_object$data
  out    <- psfm_object$out
  form   <- Formula::as.Formula(psfm_object$formula)
  n_rhs  <- length(form)[2]   ## number of RHS parts (1 = x only, 2 = x|z, 3 = x|z|h)
  
  if (!(individual %in% names(data))) {
    stop("Column '", individual, "' (the `individual` argument) was not found in psfm_object$data.")
  }
  
  ## ---- 1. Resolve h_type if "auto" -----------------------------------------
  if (h_type == "auto") {
    h_type <- if (n_rhs >= 3) "parametric" else "none"
  }
  if (h_type == "parametric" && n_rhs < 3) {
    stop("h_type = 'parametric' requires a 3-part formula (y ~ x | z | h), ",
         "but the model formula only has ", n_rhs, " RHS part(s).")
  }
  
  ## ---- 2. Build model matrices for each formula block ----------------------
  ## model.matrix() on each Formula part gives us the exact covariate columns
  ## (including/excluding intercept) the model itself used -- no need to hand
  ## list column names.
  data_x <- model.matrix(form, data = data, rhs = 1)
  data_z <- model.matrix(form, data = data, rhs = 2)
  data_h <- if (n_rhs >= 3) model.matrix(form, data = data, rhs = 3) else NULL
  
  Kx <- ncol(data_x)
  Kz <- ncol(data_z)
  Kh <- if (h_type == "parametric") ncol(data_h) else 0L
  
  ## ---- 3. Slice parameter rows out of $out ----------------------------------
  ## Row layout: sigv, sigr, [x-block], [z-block], [h-block (if parametric)]
  n_par <- nrow(out)
  expected_n_par <- 2 + Kx + Kz + Kh
  if (h_type == "scalar") expected_n_par <- expected_n_par + 1L  ## one extra row for sigma_h
  
  if (n_par != expected_n_par) {
    stop("Row count of psfm_object$out (", n_par, ") does not match the expected layout ",
         "(2 + Kx + Kz", if (h_type != "none") " + Kh" else "", " = ", expected_n_par, "). ",
         "Check that `h_type` matches how this model was actually specified.")
  }
  
  sigv_row  <- 1
  sigr_row  <- 2
  x_rows    <- (2 + 1):(2 + Kx)
  z_rows    <- (2 + Kx + 1):(2 + Kx + Kz)
  h_rows    <- if (h_type == "parametric") {
    (2 + Kx + Kz + 1):(2 + Kx + Kz + Kh)
  } else if (h_type == "scalar") {
    2 + Kx + Kz + 1
  } else {
    integer(0)
  }
  
  beta_x_hat <- out[x_rows, 1]
  beta_z_hat <- out[z_rows, 1]
  beta_h_hat <- if (h_type != "none") out[h_rows, 1] else NULL
  sigv_hat   <- out[sigv_row, 1]
  sigr_hat   <- out[sigr_row, 1]
  
  ## ---- 4. Panel bookkeeping: number of time periods per individual ---------
  ids        <- data[[individual]]
  uniq_ids   <- unique(ids)
  n_id       <- length(uniq_ids)
  timez      <- as.integer(table(factor(ids, levels = uniq_ids)))  ## preserves uniq_ids order
  
  n_obs      <- nrow(data)
  
  ## Response variable name. NOTE: all.vars() is not S3-generic, so calling it
  ## directly on a `Formula` object does not dispatch correctly and silently
  ## returns nothing usable. We instead ask Formula for just the LHS as a
  ## plain base-R formula (lhs = 1, rhs = 0), which all.vars() *does* handle.
  y_name <- all.vars(formula(form, lhs = 1, rhs = 0))[1]
  
  ## ---- 5. Set up output containers ------------------------------------------
  ## NOTE: $H is the time-invariant, individual-specific inefficiency score --
  ## one value PER INDIVIDUAL (length = n_id), NOT one value per observation
  ## (n_obs) like $U. boot_eff_h must therefore be sized by n_id, not n_obs.
  ## We derive its width from length(psfm_object$H) directly (rather than just
  ## trusting our own n_id count) so a genuine mismatch is caught immediately
  ## with a clear error instead of silently misaligning bootstrap draws.
  ## This check is unconditional on h_type: $H/boot_eff_h are always collected,
  ## even when h_type = "none" (a model can report individual effects without
  ## treating them as a parametric/scalar `h` term in the bootstrap DGP).
  n_h <- length(psfm_object$H)
  if (n_h != n_id) {
    stop("length(psfm_object$H) (", n_h, ") does not match the number of unique ",
         "individuals implied by the `individual` column (", n_id, "). ",
         "Check that `individual` is correct and that $H is one value per individual.")
  }

  boot_par             <- matrix(0, nrow = BOOT, ncol = n_par + 2)
  colnames(boot_par)   <- c(rownames(out), "loglik", "hours")
  boot_eff             <- matrix(0, nrow = BOOT, ncol = n_obs)
  boot_eff_h           <- matrix(0, nrow = BOOT, ncol = n_h)
  rownames(boot_par)   <- rownames(boot_eff) <- rownames(boot_eff_h) <- seq_len(BOOT)
  colnames(boot_eff_h) <- uniq_ids
  
  ## ---- 6. The per-replication worker function -------------------------------
  boot_one <- function(b) {
    
    set.seed(b + seed_offset)
    
    ## v: classic symmetric noise
    v <- rnorm(n_obs, 0, sigv_hat)
    
    ## u: half-normal with covariate-driven sigma (z-block), as in the original
    sigma_u <- sqrt(exp(as.vector(data_z %*% beta_z_hat)))
    u       <- abs(rnorm(n_obs, 0, sigma_u))
    
    ## r: individual random effect, constant within individual, repeated over time
    r_i <- rnorm(n_id, 0, sigr_hat)
    r   <- rep(r_i, times = timez)
    
    ## h: individual-specific inefficiency, constant within individual
    h <- switch(h_type,
                "none" = 0,
                "scalar" = {
                  sigma_h <- sqrt(exp( unname(beta_h_hat) ))  ## single scalar sigma_h, same for all individuals
                  h_i <- abs(rnorm(n_id, 0, sigma_h))
                  rep(h_i, times = timez)
                },
                "parametric" = {
                  ## one row of data_h per individual is needed; assume data_h is constant
                  ## within individual (h is a time-invariant effect) -- take the first
                  ## row observed for each individual to build the per-individual sigma.
                  first_idx  <- match(uniq_ids, ids)
                  sigma_h_i  <- sqrt(exp(as.vector(data_h[first_idx, , drop = FALSE] %*% beta_h_hat)) )
                  h_i        <- abs(rnorm(n_id, 0, sigma_h_i))
                  rep(h_i, times = timez)
                }
    )
    
    ## ---- simulate new response ----
    ## Sign convention matches sfa::psfm()'s `inefdec` argument:
    ##   inefdec = FALSE -> cost function:       inefficiency RAISES cost  (+u, +h)
    ##   inefdec = TRUE  -> production function: inefficiency LOWERS output (-u, -h)
    data_b <- data
    if(inefdec == FALSE){
      data_b[[y_name]] <- as.vector(data_x %*% beta_x_hat) + v + u + r + h  
    }else{
      data_b[[y_name]] <- as.vector(data_x %*% beta_x_hat) + v - u + r - h
    }
    
    ## ---- re-estimate ----
    MOD <- tryCatch(
      sfa::psfm(
        formula       = form,
        model_name    = psfm_object$model_name,
        data          = data_b,
        maxit.psoptim = maxit.psoptim,
        maxit.bobyqa  = maxit.bobyqa,
        maxit.optim   = maxit.optim, 
        individual    = individual, 
        inefdec       = inefdec,
        PSopt         = TRUE,
        optHessian    = FALSE,
        rand.gtre     = rand.gtre,
        rand.psoptim  = rand.psoptim, 
        start_val     = out[, 1]
      ),
      error = function(e) e
    )
    
    if (inherits(MOD, "error")) {
      return(list(sim = b, ok = FALSE, msg = conditionMessage(MOD)))
    }

    ## Defensive shape check: convert a malformed $U/$H from this particular
    ## re-estimation into a per-replication failure (recorded and skipped)
    ## rather than a fatal error that aborts the entire bootstrap after all
    ## other replications have already run.
    if (length(MOD$U) != n_obs) {
      return(list(sim = b, ok = FALSE,
                  msg = paste0("length(MOD$U) = ", length(MOD$U),
                               " does not match n_obs = ", n_obs)))
    }
    if (length(MOD$H) != n_h) {
      return(list(sim = b, ok = FALSE,
                  msg = paste0("length(MOD$H) = ", length(MOD$H),
                               " does not match n_h = ", n_h)))
    }

    list(
      sim      = b,
      ok       = TRUE,
      eff      = MOD$U,
      eff_h    = MOD$H,
      par_row  = c(MOD$out[, 1], MOD$opt$value, as.numeric(MOD$total_time, units = "hours"))
    )
  }
  
  ## ---- 7. Run in parallel ----------------------------------------------------
  cl <- parallel::makeCluster(numCores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c("data", "data_x", "data_z", "data_h", "form", "psfm_object", "out",
                "beta_x_hat", "beta_z_hat", "beta_h_hat", "sigv_hat", "sigr_hat",
                "h_type", "h_rows", "timez", "uniq_ids", "ids", "n_id", "n_obs", "n_h",
                "y_name", "individual", "maxit.psoptim", "seed_offset"),
    envir = environment()
  )
  for (p in pkgs) {
    parallel::clusterCall(cl, function(pkg) library(pkg, character.only = TRUE), p)
  }
  
  if (requireNamespace("pbapply", quietly = TRUE)) {
    results <- pbapply::pblapply(X = seq_len(BOOT), FUN = boot_one, cl = cl)
  } else {
    message("Package 'pbapply' not installed -- running without a progress bar. ",
            "Install it (install.packages(\"pbapply\")) to see live progress.")
    results <- parallel::parLapply(cl = cl, X = seq_len(BOOT), fun = boot_one)
  }
  
  ## ---- 8. Assemble results, flag failures -----------------------------------
  failures <- vapply(results, function(r) !isTRUE(r$ok), logical(1))
  
  for (b in seq_len(BOOT)) {
    if (!failures[b]) {
      boot_par[b, ]   <- results[[b]]$par_row
      boot_eff[b, ]   <- results[[b]]$eff
      boot_eff_h[b, ] <- results[[b]]$eff_h
    } else {
      boot_par[b, ]   <- NA
      boot_eff[b, ]   <- NA
      boot_eff_h[b, ] <- NA
    }
  }
  
  failed_idx <- which(failures)
  if (length(failed_idx) > 0) {
    failed_msgs <- vapply(results[failed_idx], function(r) r$msg, character(1))
    warning(length(failed_idx), " of ", BOOT, " bootstrap replications failed to ",
            "re-estimate and were set to NA.\n",
            paste0("  [b=", failed_idx, "]: ", failed_msgs, collapse = "\n"))
  }
  
  ## ---- 9. Bootstrap SEs / t-values for each parameter in $out ---------------
  boot_se   <- apply(boot_par[, seq_len(n_par), drop = FALSE], 2, sd, na.rm = TRUE)
  boot_tval <- out[, 1] / boot_se
  names(boot_se) <- names(boot_tval) <- rownames(out)
  
  result <- list(
    boot_par   = boot_par,
    boot_eff   = boot_eff,
    boot_eff_h = boot_eff_h,
    se         = boot_se,
    tval       = boot_tval,
    failures   = if (length(failed_idx) > 0) failed_idx else NULL
  )
  
  if (write_back) {
    model_out <- psfm_object
    model_out$out[, 2] <- boot_se
    model_out$out[, 3] <- boot_tval
    result$model <- model_out
  }
  
  result
}
