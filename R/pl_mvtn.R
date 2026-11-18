## Pitt and Lee (1981, JoE) Model III: the multivariate truncated normal panel
## likelihood from their Appendix 2.

## Pieces of the equicorrelated system. Returns NULL on an inadmissible
## parameter vector, which every caller treats as "return the penalty".
.pl_mvtn_parts <- function(sig_u, rho, sig_v, BigT) {
  if (any(!is.finite(c(sig_u, rho, sig_v))) || sig_u <= 0 || sig_v <= 0) {
    return(NULL)
  }
  ## rho must keep Sigma positive definite: -1/(T-1) < rho < 1
  if (rho <= -1 / (BigT - 1) + 1e-8 || rho >= 1 - 1e-8) {
    return(NULL)
  }

  a <- sig_u^2 * (1 - rho) ## Sigma = a I + b J
  b <- sig_u^2 * rho
  if (a <= 0 || a + BigT * b <= 0) {
    return(NULL)
  }
  logdet_Sig <- (BigT - 1) * log(a) + log(a + BigT * b)

  ## Sigma^{-1} = (1/a) I - cc J
  cc <- b / (a * (a + BigT * b))

  ## Q^{-1} = Sigma^{-1} + I/sigma_v^2 = a2 I - cc J
  a2 <- 1 / a + 1 / sig_v^2
  den <- a2 - BigT * cc
  if (!is.finite(a2) || a2 <= 0 || !is.finite(den) || den <= 0) {
    return(NULL)
  }
  logdet_Q <- -((BigT - 1) * log(a2) + log(den))

  ## Q = (1/a2) I + qq J
  qq <- cc / (a2 * den)
  Q <- diag(1 / a2, BigT) + qq
  if (any(!is.finite(Q))) {
    return(NULL)
  }

  list(Q = Q, logdet_Sig = logdet_Sig, logdet_Q = logdet_Q,
       a = a, b = b, Sigma = diag(a, BigT) + b, r_Q = 1 / a2, q_Q = qq)
}

## Orthant probability log Pr(w <= 0) for w ~ N(m, r I + q J), one per ROW of
## M (A100). Every covariance this likelihood meets is equicorrelated.
## q >= 0: w_t = m_t + sqrt(r) z_t + sqrt(q) c with c ~ N(0, 1), so
##   Pr = INT phi(c) prod_t Phi(-(m_t + sqrt(q) c) / sqrt(r)) dc,
## whose log integrand is concave: adaptive Gauss-Hermite at the mode (Newton,
## vectorised over rows), scaled by the curvature there. Exact to ~1e-12
## relative, including far in the tail where the old mnormt::sadmvn() at its
## defaults ("accuracy non achieved", status discarded) was off by ~3e-3.
## q < 0 has no real one-factor form: mvtnorm's Miwa rule (deterministic) up
## to T = 20, GenzBretz at a tight tolerance beyond.
.equi_log_orthant <- function(M, r, q, n_nodes = 48L) {
  if (!is.matrix(M)) M <- matrix(M, nrow = 1L)
  n <- nrow(M); Tn <- ncol(M)
  if (!(is.finite(r) && r > 0 && is.finite(q) && r + Tn * q > 0)) {
    return(rep(NA_real_, n))
  }
  if (q < 0) {
    S <- diag(r, Tn) + q
    alg <- if (Tn <= 20L) mvtnorm::Miwa(steps = 512L) else
      mvtnorm::GenzBretz(maxpts = 1e6, abseps = 1e-12, releps = 1e-8)
    .rng_state <- .rng_snapshot(); on.exit(.rng_restore(.rng_state), add = TRUE)
    set.seed(20261118L)
    return(apply(M, 1L, function(m) {
      p <- tryCatch(as.numeric(mvtnorm::pmvnorm(upper = rep(0, Tn), mean = m,
        sigma = S, algorithm = alg)), error = function(e) NA_real_)
      if (is.finite(p) && p > 0) log(p) else NA_real_
    }))
  }
  sr <- sqrt(r)
  A <- -M / sr
  b <- -sqrt(q) / sr
  g <- function(cc) -cc^2 / 2 + rowSums(pnorm(A + b * cc, log.p = TRUE))
  cc <- rep(0, n); gc <- g(cc)
  for (it in 1:100) {
    X <- A + b * cc
    Mi <- exp(dnorm(X, log = TRUE) - pnorm(X, log.p = TRUE))
    g1 <- -cc + b * rowSums(Mi)
    g2 <- -1 - b^2 * rowSums(Mi * (X + Mi))
    step <- -g1 / g2
    for (h in 1:40) {                  # damped: never step downhill
      cn <- cc + step; gn <- g(cn)
      bad <- !(gn >= gc - 1e-12 * abs(gc))
      if (!any(bad)) break
      step[bad] <- step[bad] / 2
    }
    cc <- cn; gc <- gn
    if (max(abs(step)) < 1e-12) break
  }
  X <- A + b * cc
  Mi <- exp(dnorm(X, log = TRUE) - pnorm(X, log.p = TRUE))
  sd_c <- 1 / sqrt(1 + b^2 * rowSums(Mi * (X + Mi)))
  gh <- .gauss_hermite_nodes(n_nodes)
  L <- vapply(seq_along(gh$nodes), function(k)
    g(cc + sqrt(2) * sd_c * gh$nodes[k]) + gh$nodes[k]^2 + log(gh$weights[k]), numeric(n))
  if (!is.matrix(L)) L <- matrix(L, nrow = n)
  mx <- apply(L, 1L, max)
  mx + log(rowSums(exp(L - mx))) + log(sqrt(2) * sd_c) - 0.5 * log(2 * pi)
}

## Posterior scores for the firm's w ~ N(mu_i, r I + q J) truncated to w <= 0
## (inefficiency = -w), one row per firm (A100):
##   E[exp(w_t)]  = exp(mu_t + (r+q)/2) Pr(mu + Q e_t) / Pr(mu)   (truncated MGF)
##   E[w_t]       = mu_t - sum_s Q_ts F_s / Pr(mu)                (Tallis 1961)
##   F_s = phi(0; mu_s, r+q) Pr_{T-1}(w_-s <= 0 | w_s = 0), whose conditional
##   law is again equicorrelated: mean mu_-s - q/(r+q) mu_s, cov r I + qr/(r+q) J.
.pl_mvtn_scores <- function(MU, r, q) {
  Tn <- ncol(MU); v <- r + q
  lP <- .equi_log_orthant(MU, r, q)
  BC <- U <- matrix(NA_real_, nrow(MU), Tn)
  qc <- q * r / v
  LF <- matrix(NA_real_, nrow(MU), Tn)
  for (t in seq_len(Tn)) {
    sh <- MU + q; sh[, t] <- sh[, t] + r          # mu + Q e_t
    BC[, t] <- exp(MU[, t] + v / 2 + .equi_log_orthant(sh, r, q) - lP)
    Mc <- MU[, -t, drop = FALSE] - (q / v) * MU[, t]
    LF[, t] <- dnorm(0, MU[, t], sqrt(v), log = TRUE) + .equi_log_orthant(Mc, r, qc) - lP
  }
  FP <- exp(LF)                                   # F_s / Pr
  Ew <- MU - (r * FP + q * rowSums(FP))
  list(u_hat = pmax(-Ew, 0), exp_u_hat = pmin(pmax(BC, 0), 1))
}

## The out-of-domain penalty. NOT .SFA_CONSTANTS$MAX_VALUE, which is
## .Machine$double.xmax: this is the one likelihood in the package handed
## straight to stats::optim() with no derivative-free stage in front of it, so
## L-BFGS-B differences the penalty against an ordinary value and
## (1.8e308 - 239)/1e-8 overflows to Inf -- "non-finite value supplied by
## optim" before the fit starts. ^0.1 is the same damping psfm() already
## applies where an optimizer sees the penalty directly.
.PL_MVTN_PEN <- function() .SFA_CONSTANTS$MAX_VALUE^0.1

## Negative summed log-likelihood. `par` is (sig_v, sig_u, rho, beta...);
## eps_list holds one T-vector of y - X beta per firm, rebuilt each call.
## With per_obs = TRUE it returns the POSITIVE per-FIRM contributions instead.
.pl_mvtn_nll <- function(par, Y_list, X_list, BigT, inefdec_n = 1,
                         per_obs = FALSE) {
  sig_v <- par[1]
  sig_u <- par[2]
  rho <- par[3]
  beta <- par[-(1:3)]

  nfirm <- length(Y_list)
  ## The penalty must keep its length when per-firm contributions are asked
  ## for: estfun() central-differences this, and a scalar returned beside a
  ## length-nfirm vector recycles into a meaningless score column.
  .bail <- function() {
    if (isTRUE(per_obs)) {
      rep(-.PL_MVTN_PEN() / nfirm, nfirm)
    } else {
      .PL_MVTN_PEN()
    }
  }

  pr <- .pl_mvtn_parts(sig_u, rho, sig_v, BigT)
  if (is.null(pr)) {
    return(.bail())
  }
  logP0 <- .equi_log_orthant(rep(0, BigT), pr$a, pr$b)
  if (!is.finite(logP0)) {
    return(.bail())
  }

  Qi <- pr$Q
  s2 <- sig_v^2
  const <- -(BigT / 2) * log(2 * pi) - BigT * log(sig_v) -
    0.5 * pr$logdet_Sig + 0.5 * pr$logdet_Q - logP0

  ## u <= 0 is PL's sign convention, i.e. inefficiency SUBTRACTS from a
  ## production frontier. inefdec_n flips it for a cost frontier.
  E <- t(vapply(seq_len(nfirm), function(i)
    inefdec_n * (Y_list[[i]] - as.numeric(X_list[[i]] %*% beta)), numeric(BigT)))
  QE <- E %*% Qi                                  # Q symmetric
  MU <- QE / s2
  quad <- rowSums(E^2) / s2 - rowSums(E * QE) / s2^2
  lP <- .equi_log_orthant(MU, pr$r_Q, pr$q_Q)
  ll_i <- const - 0.5 * quad + lP
  if (any(!is.finite(ll_i))) {
    if (!isTRUE(per_obs)) return(.PL_MVTN_PEN())
    ll_i[!is.finite(ll_i)] <- -.PL_MVTN_PEN() / nfirm
  }
  if (isTRUE(per_obs)) {
    ## Positive, one entry per FIRM, non-finite floored rather than dropped.
    ll_i[!is.finite(ll_i)] <- -.PL_MVTN_PEN() / nfirm
    return(ll_i)
  }
  tot <- sum(ll_i)
  if (!is.finite(tot)) {
    return(.PL_MVTN_PEN())
  }
  -tot
}


## Self-contained fitting routine, called from psfm() before the shared
## start_panel() machinery.
.psfm_pl_mvtn <- function(data, y_var, x_vars_vec, individual, inefdec_n,
                          maxit.optim, Method, optHessian, start_val,
                          verbose, call, formula, keep_objective = FALSE) {
  id_chr <- as.character(data[, individual])
  indiv <- unique(id_chr)
  gid <- match(id_chr, indiv)
  t_i <- as.numeric(table(gid))
  if (length(unique(t_i)) != 1L) {
    stop("psfm(model_name = \"PL80_MVTN\") requires a BALANCED panel: Sigma is ",
      "a single T x T matrix shared by every firm. Found T ranging from ",
      min(t_i), " to ", max(t_i), ". Use model_name = \"PL80\" for the ",
      "time-invariant model, which handles unbalanced panels.",
      call. = FALSE
    )
  }
  BigT <- t_i[1]
  if (BigT < 2L) {
    stop("psfm(model_name = \"PL80_MVTN\") needs at least two periods per firm; ",
      "with T = 1 there is no cross-period dependence for Sigma to describe.",
      call. = FALSE
    )
  }

  ## "(Intercept)" is a NAME in x_vars_vec, not a column of `data` -- the same
  ## thing the GTRE_FML branch works around by setdiff()-ing it out.
  xn <- x_vars_vec
  xcols <- setdiff(xn, "(Intercept)")
  Xall <- as.matrix(as.data.frame(data)[, xcols, drop = FALSE])
  if ("(Intercept)" %in% xn) {
    Xall <- cbind(`(Intercept)` = 1, Xall)
    xn <- c("(Intercept)", xcols)
  } else {
    xn <- xcols
  }
  Yall <- as.numeric(as.data.frame(data)[, y_var])
  Y_list <- split(Yall, gid)
  X_list <- lapply(split(seq_along(gid), gid), function(ix) Xall[ix, , drop = FALSE])

  ## Start values: pooled OLS for beta, and the usual moment split of the
  ## residual variance between the two components.
  ols <- stats::lm.fit(Xall, Yall)
  b0 <- as.numeric(ols$coefficients)
  b0[!is.finite(b0)] <- 0
  r <- as.numeric(ols$residuals)
  s2 <- stats::var(r)
  start_v <- if (isTRUE(is.numeric(start_val)) && length(start_val) == length(b0) + 3L) {
    as.numeric(start_val)
  } else {
    c(sqrt(max(0.5 * s2, 1e-3)), sqrt(max(0.5 * s2, 1e-3)), 0.3, b0)
  }

  Start.Time <- start.time()
  np <- length(start_v)
  lower <- c(1e-3, 1e-3, -1 / (BigT - 1) + 0.02, rep(-Inf, np - 3))
  upper <- c(Inf, Inf, 0.95, rep(Inf, np - 3))

  ## Three stages, not one. The orthant probability inside the likelihood is
  ## computed to limited precision, so its finite-difference gradient is
  ## noisy and L-BFGS-B alone stops at points that are not maxima -- while
  ## still reporting convergence = 0. A derivative-free stage between two
  ## L-BFGS-B runs is the same scaffold every other psfm() model uses. It
  ## leaves already-converged fits untouched and recovers several
  ## log-likelihood points where the single stage quit early.
  .nll <- function(p) {
    .pl_mvtn_nll(p,
      Y_list = Y_list, X_list = X_list, BigT = BigT, inefdec_n = inefdec_n
    )
  }
  o0 <- stats::optim(
    par = start_v, fn = .nll, method = Method,
    lower = lower, upper = upper, hessian = FALSE,
    control = list(maxit = maxit.optim, trace = if (isTRUE(verbose)) 1 else 0)
  )
  polished <- tryCatch(
    {
      bo <- minqa::bobyqa(o0$par, .nll,
        lower = lower, upper = pmin(upper, .Machine$double.xmax^0.25),
        control = list(maxfun = 5000, iprint = if (isTRUE(verbose)) 1 else 0)
      )
      if (is.finite(bo$fval) && bo$fval <= o0$value) bo$par else o0$par
    },
    error = function(e) o0$par
  )
  opt <- stats::optim(
    par = polished, fn = .nll, method = Method,
    lower = lower, upper = upper, hessian = isTRUE(optHessian),
    control = list(maxit = maxit.optim, trace = if (isTRUE(verbose)) 1 else 0)
  )
  ## Never return a worse point than the first stage found.
  if (!is.finite(opt$value) || opt$value > o0$value) {
    opt <- stats::optim(
      par = o0$par, fn = .nll, method = Method,
      lower = lower, upper = upper, hessian = isTRUE(optHessian),
      control = list(maxit = 1L)
    )
  }
  total_time <- start.time() - Start.Time

  par_names <- c("sigv", "sigu", "rho", xn)
  out <- matrix(NA_real_, 3, np, dimnames = list(c("par", "st_err", "t-val"), par_names))
  out[1, ] <- opt$par
  se <- rep(NA_real_, np)
  if (isTRUE(optHessian) && !is.null(opt$hessian)) {
    V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) {
      se <- sqrt(diag(V))
    }
  }
  out[2, ] <- se
  out[3, ] <- out[1, ] / se

  ## The exact posterior scores (A100): E[u | eps_i] and E[exp(-u) | eps_i]
  ## under the truncated equicorrelated normal, by .pl_mvtn_scores().
  pr <- .pl_mvtn_parts(opt$par[2], opt$par[3], opt$par[1], BigT)
  u_hat <- exp_u_hat <- rep(NA_real_, length(Yall))
  if (!is.null(pr)) {
    E <- t(vapply(seq_along(Y_list), function(i)
      inefdec_n * (Y_list[[i]] - as.numeric(X_list[[i]] %*% opt$par[-(1:3)])), numeric(BigT)))
    sc <- .pl_mvtn_scores(E %*% pr$Q / opt$par[1]^2, pr$r_Q, pr$q_Q)
    for (i in seq_along(Y_list)) {
      u_hat[gid == i] <- sc$u_hat[i, ]
      exp_u_hat[gid == i] <- sc$exp_u_hat[i, ]
    }
  }

  ## t(out), not out. Every other entry point stores the p x 3 transpose, and
  ## the documented access is fit$out[, "par"] -- storing the 3 x p matrix
  ## here made that error with "subscript out of bounds" on this model alone.
  results <- list(
    out = t(out), opt = opt, total_time = total_time, start_v = start_v,
    model_name = "PL80_MVTN", formula = formula,
    coefficients = out[1, ], std.errors = out[2, ], t.values = out[3, ],
    u_hat = u_hat, exp_u_hat = exp_u_hat,
    BigT = BigT, call = call, data = data
  )
  class(results) <- "sfareg"
  results$inefdec <- isTRUE(inefdec_n > 0)
  if (isTRUE(keep_objective)) {
    ## Contributions are per FIRM, so n_units -- not nobs() -- is what
    ## bread() must scale by.
    results$objective <- function(par, per_obs = FALSE) {
      .pl_mvtn_nll(par,
        Y_list = Y_list, X_list = X_list, BigT = BigT,
        inefdec_n = inefdec_n, per_obs = per_obs
      )
    }
    results$n_units <- length(Y_list)
  }
  results
}
