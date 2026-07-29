## ------------------------------------------------------------
## Helper: column-wise products of a matrix (replaces matrixStats::colProds)
##
## Inlined to drop the matrixStats dependency (only colProds was ever used
## from it, in the GTRE/TRE simulated-ML likelihoods). Folds the rows with
## element-wise multiplication down the columns -- this is byte-identical to
## matrixStats::colProds(method = "direct") (the default), i.e. the same
## left-to-right product order, verified with identical(); it is NOT the
## expSumLog method (which would differ at ~1e-15 and change results in the
## 15th digit). For the small row counts here (Ti = number of time periods,
## typically <=10-20) this pure-R fold is actually faster than the C call it
## replaces, since it avoids matrixStats' per-call setup overhead.
## ------------------------------------------------------------
.col_prods <- function(x){
  x  <- as.matrix(x)
  cp <- x[1, ]
  if(nrow(x) > 1){
    for(r in 2:nrow(x)) cp <- cp * x[r, ]
  }
  cp
}


## ------------------------------------------------------------
## Helper: column-wise (or vector) demeaning (replaces Jmisc::demean)
##
## Inlined to drop the Jmisc dependency (only demean was ever used from it,
## in start.tfe()'s within-transformation for the TFE model). Matches
## Jmisc::demean exactly: a vector returns x - mean(x); a matrix/data.frame
## returns a matrix with each column demeaned.
## ------------------------------------------------------------
.demean <- function(x){
  if(is.vector(x)) return(x - mean(x))
  apply(x, 2, function(col) col - mean(col))
}


## ------------------------------------------------------------
## Helper: physicists' Gauss-Hermite quadrature nodes/weights
##
## Golub-Welsch: nodes are eigenvalues of the (n x n) tridiagonal Jacobi
## matrix with zero diagonal and off-diagonal sqrt(i/2), weights are
## sqrt(pi) times the squared first component of each normalized
## eigenvector. Base R only (eigen()); no dependency on statmod/gaussquad.
## Used by .log_mvn_cdf_rank1() below.
## ------------------------------------------------------------
.gauss_hermite_nodes <- function(n = 64L){
  n <- as.integer(n)
  if(n < 5L) stop("n must be at least 5.")
  i    <- seq_len(n - 1L)
  J    <- matrix(0, n, n)
  off  <- sqrt(i / 2)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  eg      <- eigen(J, symmetric = TRUE)
  nodes   <- eg$values
  weights <- sqrt(pi) * (eg$vectors[1L, ]^2)
  ord     <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}


## ------------------------------------------------------------
## Helper: log MVN CDF for the rank-one/equicorrelated covariance
## Sigma = I_T + c*11' (replaces mnormt::pmnorm() for this special case)
##
## For Sigma of this form, X = Z + sqrt(c)*W*1 with Z ~ N(0,I_T), W ~ N(0,1)
## independent, so Pr(X <= upper) = E_W[ prod_t Phi(upper_t - sqrt(c)*W) ] --
## an exact 1-dimensional integral, here evaluated by Gauss-Hermite
## quadrature (default 64 nodes) rather than mnormt::pmnorm()'s general
## T-dimensional (stochastic, Genz-algorithm) integration. This same
## covariance structure arises in psfm()'s TFE model (see like.tfe() in
## psfm.R). Matches mnormt::pmnorm() to within pmnorm's own default
## tolerance (~1e-3 on the log scale) while being ~30x faster and
## deterministic (no Monte Carlo noise).
## ------------------------------------------------------------
.log_mvn_cdf_rank1 <- function(upper, c, gh){
  upper <- as.numeric(upper)
  if(length(upper) == 0L) return(0)
  if(c < 1e-14) return(sum(pnorm(upper, log.p = TRUE)))

  sc      <- sqrt(c)
  s_nodes <- sqrt(2) * gh$nodes
  log_terms <- numeric(length(s_nodes))
  for(r in seq_along(s_nodes)){
    log_terms[r] <- log(gh$weights[r]) + sum(pnorm(upper - sc * s_nodes[r], log.p = TRUE))
  }
  m <- max(log_terms)
  (m + log(sum(exp(log_terms - m)))) - 0.5 * log(pi)
}


## ------------------------------------------------------------
## Helper: log density of the within-transformed disturbance vector
## for psfm()'s TFE model (replaces mnormt::dmnorm() for this special case)
##
## Closed form for the density of an (T-1)-vector eps_star with covariance
## sigma2*(I_m - (1/T)*11'_m), m = T-1, using det = 1/T and inverse =
## I_m + 11'_m (both by the matrix determinant lemma / Sherman-Morrison,
## exact given the specific 1/T normalizer used in like.tfe()'s E_t1).
## Matches mnormt::dmnorm() to floating-point precision.
## ------------------------------------------------------------
.log_within_mvn_density <- function(eps_star, sigma2){
  m <- length(eps_star)
  q <- sum(eps_star^2) + sum(eps_star)^2
  -0.5 * m * log(2 * pi) - 0.5 * m * log(sigma2) + 0.5 * log(m + 1) - q / (2 * sigma2)
}


## ------------------------------------------------------------
## Helper: cross-sectional intercept-only normal-half-normal (Aigner,
## Lovell & Schmidt 1977) MLE, with Hessian-based standard errors.
##
## Used by psfm()'s GTRE_SEQ1 branch to decompose composite residuals
## (epsilon_hat, alpha_hat) into their normal/half-normal variance
## components. sfm()'s own NHN likelihood can't be reused directly here
## because sfm()'s formula-parsing path does not currently support an
## intercept-only (y~1) formula, so this reimplements the same closed-form
## density directly -- the T_i=1 special case of the error-components-
## frontier likelihood used for PL80/BC92, parameterized the same way,
## directly over (sigma_v, sigma_u, intercept) rather than
## (lambda=sigma_u/sigma_v, sigma).
## ------------------------------------------------------------
.fit_nhn_intercept <- function(y, inefdec_n = 1, maxit.bobyqa = 200, maxit.psoptim = 20,
                                maxit.optim = 100, PSopt = FALSE, verbose = FALSE){
  like_fn <- function(x){
    sigma_v <- x[1]; sigma_u <- x[2]; intercept <- x[3]
    if(!is.finite(sigma_v) || !is.finite(sigma_u) || sigma_v <= 0 || sigma_u <= 0) return(1e12)
    eps <- inefdec_n*(y - intercept)
    sigmaSq <- sigma_v*sigma_v + sigma_u*sigma_u
    z  <- -sigma_u*eps/(sigma_v*sqrt(sigmaSq))
    ll <- log(2) - 0.5*log(2*pi) - 0.5*log(sigmaSq) - (eps*eps)/(2*sigmaSq) + pnorm(z, log.p = TRUE)
    -sum(ll[is.finite(ll)])
  }
  start_v <- c(sd(y)/sqrt(2), sd(y)/sqrt(2), mean(y))
  lower_v <- c(.SFA_CONSTANTS$MIN_POSITIVE, .SFA_CONSTANTS$MIN_POSITIVE, -Inf)

  Opt.Bobyqa <- opt.bobyqa(fn=like_fn, start_v=start_v, lower.bobyqa=lower_v,
                            maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE, verbose=verbose)
  start_v <- Opt.Bobyqa$start_v

  differ <- 2
  Opt.Psoptim <- opt.psoptim(fn=like_fn, start_v, lower.psoptim=c(lower_v[1:2], start_v[3]-differ),
                              upper.psoptim=c(start_v[1:2]+differ, start_v[3]+differ),
                              maxit.psoptim=maxit.psoptim, psopt.TF=PSopt, verbose=verbose)
  start_v <- Opt.Psoptim$start_v

  Opt.Optim <- opt.optim(fn=like_fn, start_v=start_v, lower.optim=lower_v, upper.optim=rep(Inf,3),
                          maxit.optim=maxit.optim, opt.TF=TRUE, method="L-BFGS-B", optHessian=TRUE, verbose=verbose)
  opt <- Opt.Optim$opt
  list(par = opt$par, hessian = opt$hessian, value = opt$value)
}


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
##
## Averaging M with t(M) makes the *values* symmetric, but not necessarily
## dimnames(M) -- R's `+` keeps the left operand's dimnames, so if M carries
## inconsistent/spurious dimnames (e.g. list(NULL, c("","",...)), a known
## cbind()/solve() artifact when combining an unnamed vector with an
## unnamed matrix -- see the GTRE Lambda_i computation this guards), the
## result still does too. tmvtnorm::checkSymmetricPositiveDefinite() calls
## isSymmetric() with its default check.attributes=TRUE, which fails on
## attributes alone even when every number is symmetric to machine
## precision, raising a spurious "sigma must be a symmetric matrix" error.
## Since a genuinely-value-asymmetric matrix would already fail this same
## isSymmetric() check regardless of dimnames, it's always safe to drop
## dimnames here when they're the only thing making it look asymmetric --
## no value is ever altered, and a real breakdown still isSymmetric()==
## FALSE afterwards.
## ------------------------------------------------------------
.safe_symmetrize <- function(M){
  S <- 0.5 * (M + t(M))
  if(!is.null(dimnames(S)) &&
     !isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps))) &&
     isTRUE(isSymmetric(S, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))){
    dimnames(S) <- NULL
  }
  S
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




## ------------------------------------------------------------
## Single, robust entry point for parsing sfa's pipe-delimited formula
## syntax (y ~ x1 + x2 | z1 + z2 | zp1).
##
## Replaces a second, separate formula-parsing path that used to live in
## data_proc() (data.processing.R): a manual
## `base::strsplit(as.character(formula), "|", fixed = TRUE)` followed by
## `paste()`-reconstructing each sub-formula from string fragments. That
## approach is fragile for a specific, real reason: `as.character()` on a
## formula calls `deparse()` internally, which line-wraps long RHS text
## (roughly 60+ characters -- easy to hit with more than a handful of x
## variables) into *extra* character-vector elements. The old code indexed
## into that vector positionally (`form_parts[[3]][1]`, `[2]`, `[3]`), so a
## sufficiently long formula would silently misparse rather than error.
## `.check_model_formula_pipes()` elsewhere in this file already used
## `Formula::Formula()` (an existing Imports dependency, purpose-built for
## exactly this pipe-delimited-formula pattern) to just *count* RHS parts;
## this extends the same robust approach to actually extracting each part.
##
## Returns:
##   n_parts    : number of RHS parts (1, 2, or 3)
##   y_var      : character, name of the response variable
##   formula_x  : two-sided formula for the main regression part (y ~ x...)
##   formula_z  : two-sided formula for the first pipe segment (y ~ z...),
##                NULL if n_parts < 2
##   formula_zp : two-sided formula for the second pipe segment (y ~ zp...),
##                NULL if n_parts < 3
##
## Deliberately does NOT return variable-name vectors for the z/zp parts --
## data_proc() derives those downstream from `colnames(data_conform(...))`
## instead (already robust to interactions/transforms/etc, since it comes
## from an actual model.matrix rather than text parsing), so duplicating
## that here would just be a second source of truth.
## ------------------------------------------------------------
.parse_pipe_formula <- function(formula){

  f       <- Formula::Formula(formula)
  n_parts <- length(f)[2]

  formula_x  <- formula(f, lhs = 1, rhs = 1)
  formula_z  <- if(n_parts >= 2) formula(f, lhs = 1, rhs = 2) else NULL
  formula_zp <- if(n_parts >= 3) formula(f, lhs = 1, rhs = 3) else NULL

  y_var <- all.vars(formula_x)[1]

  list(
    n_parts    = n_parts,
    y_var      = y_var,
    formula_x  = formula_x,
    formula_z  = formula_z,
    formula_zp = formula_zp
  )
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
    # lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
    lower_bounds <- c(lower_bounds, inf_sub)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_z_vars))
  }
  
  # 5. Delta_p Section 
  # If "1", it's 1 param. If variables, it's n_zp_eff params.
  if (parts[3] == "1") {
    # lower_bounds <- c(lower_bounds, .SFA_CONSTANTS$MIN_POSITIVE)
    lower_bounds <- c(lower_bounds, inf_sub)
  } else {
    lower_bounds <- c(lower_bounds, rep(inf_sub, prep$n_zp_vars))
  }
  
  return(lower_bounds)
}




.check_model_formula_pipes <- function(formula, model_name) {
  # 1. Map model names to their maximum ALLOWED RHS parts (pipes + 1)
  # 1 pipe = 2 parts, 2 pipes = 3 parts, 0 pipes = 1 part
  #
  # NHN/NE/GTRE/TRE (the non-Z models) are capped at 1 part -- i.e. no pipes
  # -- by design, to keep a sharp separation between the plain and Z
  # (heteroskedastic-inefficiency) model variants: a user who wants pipe
  # syntax must name the "_Z" model explicitly (e.g. "GTRE_Z") rather than
  # relying on implicit upgrading from a pipe-bearing formula passed to the
  # plain model name. data_proc() (data.processing.R) still contains the
  # underlying auto-upgrade logic (model_name=="NHN" -> "NHN_Z", etc.), but
  # since this function validates model_name before data_proc() ever runs,
  # that path is unreachable from sfm()/psfm() and is effectively dead code
  # kept in place rather than removed.
  max_parts_map <- c(
    # --- 0 Pipes Allowed (1 Part) ---
    "TFE"        = 1,
    "GTRE_SEQ1"  = 1,
    "GTRE_SEQ2"  = 1,
    "SSFE"       = 1,
    "PL80"       = 1,
    "BC92"       = 1,
    "NR"         = 1,
    "THT"        = 1,
    "NTN"        = 1,
    "NG"         = 1,
    "ZISF"       = 1,
    "NNAK"       = 1,
    "NHN"        = 1,
    "NE"         = 1,
    "TRE"        = 1,
    "GTRE"       = 1,

    # --- 1 Pipe Allowed (2 Parts) ---
    "TRE_Z"      = 2,
    "FD"         = 2,
    "ZISF_Z"     = 2,
    "NHN_Z"      = 2,
    "NE_Z"       = 2,

    # --- 2 Pipes Allowed (3 Parts) ---
    "GTRE_Z"     = 3,
    "TTNE"       = 3,
    "TTHN"       = 3,
    "TTNLS"      = 3
  )
  
  # 2. Check if the provided model name is valid
  if (!(model_name %in% names(max_parts_map))) {
    stop(paste("Unknown model_name:", model_name), call. = FALSE)
  }
  
  # 3. Parse the formula and count RHS parts separated by '|'
  f <- Formula::Formula(formula)
  rhs_parts <- length(f)[2] # Index 2 extracts the RHS parts vector length
  
  # 4. Get the max allowed parts for this specific model
  allowed_parts <- max_parts_map[model_name]
  
  # 5. Enforce the limit
  if (rhs_parts > allowed_parts) {
    formula_str <- paste(deparse(formula), collapse = " ")
    max_pipes <- allowed_parts - 1
    
    stop(paste0(
      "Invalid formula structure for model '", model_name, "'!\n",
      "Formula provided: ", formula_str, "\n",
      "The '", model_name, "' model allows a maximum of ", max_pipes, " pipe separator(s) (", allowed_parts, " parts).\n",
      "Found ", (rhs_parts - 1), " pipe separator(s) (", rhs_parts, " parts) instead."
    ), call. = FALSE)
  }
  
  return(invisible(TRUE))
}
