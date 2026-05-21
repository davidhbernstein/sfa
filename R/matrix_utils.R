## ------------------------------------------------------------
## Helper: construct variance-design matrices
##
## If vars is empty, return an intercept-only design so that
## the same code path nests the homoskedastic special case.
## ------------------------------------------------------------
.make_var_design <- function(data, vars, rows = NULL, int_name = "int"){
  if(length(vars) == 0){
    n_rows <- if(is.null(rows)) nrow(data) else length(rows)
    out <- matrix(1, nrow = n_rows, ncol = 1)
    colnames(out) <- int_name
    return(out)
  }
  
  out <- as.matrix(data[, vars, drop = FALSE])
  
  if(!is.null(rows)){
    out <- out[rows, , drop = FALSE]
  }
  
  out
}


## ------------------------------------------------------------
## Helper: numerical symmetry enforcement
## ------------------------------------------------------------
.safe_symmetrize <- function(M){
  0.5 * (M + t(M))
}


## ------------------------------------------------------------
## Helper: safe inverse for covariance-type matrices
##
## Strategy:
##   1. symmetrize
##   2. try Cholesky
##   3. progressively ridge if needed
##   4. fall back to solve()
## ------------------------------------------------------------
.safe_inverse <- function(M,
                          base_ridge = 1e-10,
                          ridge_mult = 10,
                          max_tries = 8,
                          name = "matrix") {
  
  M <- .safe_symmetrize(M)
  p <- nrow(M)
  I_p <- diag(p)
  
  ## First try: Cholesky path
  for(k in 0:max_tries){
    ridge <- if(k == 0) 0 else base_ridge * ridge_mult^(k - 1)
    M_try <- if(ridge == 0) M else M + ridge * I_p
    
    chol_obj <- tryCatch(chol(M_try), error = function(e) NULL)
    
    if(!is.null(chol_obj)){
      M_inv <- chol2inv(chol_obj)
      return(list(
        value   = .safe_symmetrize(M_inv),
        ridge   = ridge,
        success = TRUE,
        method  = if(ridge == 0) "chol" else "chol_ridge"
      ))
    }
  }
  
  ## Second try: direct solve fallback
  for(k in 0:max_tries){
    ridge <- if(k == 0) 0 else base_ridge * ridge_mult^(k - 1)
    M_try <- if(ridge == 0) M else M + ridge * I_p
    
    sol <- tryCatch(solve(M_try), error = function(e) NULL)
    
    if(!is.null(sol)){
      return(list(
        value   = .safe_symmetrize(sol),
        ridge   = ridge,
        success = TRUE,
        method  = if(ridge == 0) "solve" else "solve_ridge"
      ))
    }
  }
  
  stop(sprintf("Unable to invert %s even after ridging.", name))
}


## ------------------------------------------------------------
## Helper: compute Lambda_i and ARR_i robustly
##
## Lambda_i = [VEE_i^{-1} + A_i' SIG_i^{-1} A_i]^{-1}
## ARR_i    = Lambda_i A_i' SIG_i^{-1}
## ------------------------------------------------------------
.safe_linear_combo <- function(VEE_i, A_i, SIG_i,
                               base_ridge = 1e-10,
                               ridge_mult = 10,
                               max_tries = 8,
                               name = "posterior system") {
  
  invVEE <- .safe_inverse(VEE_i,
                          base_ridge = base_ridge,
                          ridge_mult = ridge_mult,
                          max_tries  = max_tries,
                          name       = "VEE_i")
  
  invSIG <- .safe_inverse(SIG_i,
                          base_ridge = base_ridge,
                          ridge_mult = ridge_mult,
                          max_tries  = max_tries,
                          name       = "SIG_i")
  
  K <- .safe_symmetrize(invVEE$value + t(A_i) %*% invSIG$value %*% A_i)
  
  invK <- .safe_inverse(K,
                        base_ridge = base_ridge,
                        ridge_mult = ridge_mult,
                        max_tries  = max_tries,
                        name       = name)
  
  ARR <- invK$value %*% t(A_i) %*% invSIG$value
  
  list(
    LAM = .safe_symmetrize(invK$value),
    ARR = ARR,
    invVEE_method = invVEE$method,
    invSIG_method = invSIG$method,
    invK_method   = invK$method,
    ridge_VEE     = invVEE$ridge,
    ridge_SIG     = invSIG$ridge,
    ridge_K       = invK$ridge
  )
}


## Function for extending user formulas
.format_formula <- function(input_val) {
  # Convert formula object to a single string
  input_string <- paste(format(input_val), collapse = "")
  
  # Split by pipes
  parts <- trimws(strsplit(input_string, "\\|")[[1]])
  
  # Replace empty parts with "1"
  parts[parts == ""] <- "1"
  
  # Pad with "1" until there are 3 components
  while (length(parts) < 3) {
    parts <- c(parts, "1")
  }
  
  # Return as a formula
 return(as.formula(paste(parts, collapse = " | ")))
}




## function to handle gtre lower bounds for sigmas
.generate_sfa_bounds <- function(input_formula, prep, inf_sub = -Inf) {
  
  # 1. Split formula into parts
  f_char <- paste(format(input_formula), collapse = "")
  parts <- trimws(strsplit(f_char, "\\|")[[1]])
  while(length(parts) < 3) parts <- c(parts, "1")
  parts[parts == ""] <- "1"
  
  # 2. Start with Variance Components (2 params)
  lower_bounds <- rep(.SFA_CONSTANTS$MIN_POSITIVE, 2)
  
  # 3. Beta Section (prep$n_x_vars params)
  lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_x_vars))
  
  # 4. Delta Section 
  # If "1", it's 1 param. If variables, it's n_z_eff params.
  if (parts[2] == "1") {
    lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_z_vars))
  }
  
  # 5. Delta_p Section 
  # If "1", it's 1 param. If variables, it's n_zp_eff params.
  if (parts[3] == "1") {
    lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_zp_vars))
  }
  
  return(lower_bounds)
}

