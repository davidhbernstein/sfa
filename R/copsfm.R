## Copula stochastic frontier: dependence between the noise and the
## inefficiency. Gap K5. See notes/code_history/copsfm.md.
##
## Every other model in this package assumes v and u are independent. That
## assumption is ubiquitous in the frontier literature and, as Smith (2008)
## argued, mostly an artefact of how easily the likelihood is then written
## down: if a farmer misjudges a seasonal shock and mis-allocates inputs in
## response, the shock and the inefficiency are the same event seen twice.
##
## The package's own taxonomy already states the general form (JSS paper, Eq 6):
##
##   f_{V,U}(v, u) = f_V(v) f_U(u) c(F_V(v), F_U(u); rho)
##
## with every existing specification setting c = 1. This file relaxes that, and
## nothing else: the marginals are the ordinary normal and half-normal.
##
## Composing, with eps = v - S*u so that v = eps + S*u,
##
##   f_eps(eps) = int_0^inf f_V(eps + S u) f_U(u) c(F_V(eps + S u), F_U(u)) du
##
## which is a one-dimensional integral per observation, done by Gauss-Legendre
## on a transformed range rather than by simulation.

## Copula densities on the unit square, in LOGS.
##
## Only families whose density can be written down exactly and checked are
## included. Each is verified in the tests two ways: it integrates to 1 over
## [0,1]^2, and at the independence parameter it returns exactly 1.
.cop_logc <- function(w1, w2, par, family) {
  switch(family,
    independent = rep(0, length(w1)),
    gaussian = {
      ## rho in (-1, 1). a = Phi^-1(w1), b = Phi^-1(w2):
      ##   c = (1-rho^2)^-1/2 exp{ [2 rho a b - rho^2 (a^2 + b^2)] / [2(1-rho^2)] }
      rho <- par
      if (abs(rho) >= 1) return(rep(NA_real_, length(w1)))
      a <- stats::qnorm(w1)
      b <- stats::qnorm(w2)
      om <- 1 - rho^2
      -0.5 * log(om) + (2 * rho * a * b - rho^2 * (a^2 + b^2)) / (2 * om)
    },
    fgm = {
      ## Farlie-Gumbel-Morgenstern: c = 1 + theta (1-2 w1)(1-2 w2), |theta| <= 1.
      ## Bounded dependence (Spearman rho is theta/3), so it cannot represent
      ## strong association -- which is exactly why it is a useful contrast to
      ## the Gaussian rather than a redundant second option.
      th <- par
      if (abs(th) > 1) return(rep(NA_real_, length(w1)))
      v <- 1 + th * (1 - 2 * w1) * (1 - 2 * w2)
      ifelse(v > 0, log(v), NA_real_)
    },
    ## NOTE: Frank, Clayton and Gumbel are deliberately ABSENT. Each is a
    ## one-line density in a reference, and a one-line density transcribed
    ## without a source to check it against is exactly how a wrong likelihood
    ## ships looking right. Adding one means: write it, verify it integrates to
    ## 1 over the unit square, and verify it returns 1 at independence -- the
    ## two checks the families below pass in test-copsfm.R.
    NA_real_
  )
}

## Independence value of each family's parameter, and its admissible range.
.cop_spec <- function(family) {
  switch(family,
    gaussian = list(par0 = 0, lo = -0.95, hi = 0.95, name = "rho"),
    fgm      = list(par0 = 0, lo = -0.99, hi = 0.99, name = "theta"),
    list(par0 = 0, lo = -0.95, hi = 0.95, name = "par")
  )
}

copsfm <- function(formula,
                   data,
                   copula = c("gaussian", "fgm"),
                   inefdec = TRUE,
                   n_nodes = 128,
                   maxit.bobyqa = 10000,
                   maxit.psoptim = 1000,
                   maxit.optim = 1000,
                   start_val = FALSE,
                   PSopt = FALSE,
                   optHessian = TRUE,
                   Method = "L-BFGS-B",
                   verbose = FALSE,
                   rand.psoptim = NULL) {
  call <- match.call()
  copula <- match.arg(copula)
  Start.Time <- Sys.time()
  cz <- .SFA_CONSTANTS

  if (missing(formula) || !inherits(formula, "formula") || length(formula) != 3L) {
    stop("copsfm(): `formula` must be a two-sided formula, e.g. y ~ x1 + x2.",
      call. = FALSE
    )
  }
  if (any(grepl("|", deparse(formula), fixed = TRUE))) {
    stop("copsfm(): `formula` must not contain a `|` segment. Determinants of ",
      "the variance are not yet supported alongside a copula.",
      call. = FALSE
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("copsfm(): `data` must be a data.frame.", call. = FALSE)
  }
  if (length(inefdec) != 1L || !is.logical(inefdec) || is.na(inefdec)) {
    stop("copsfm(): `inefdec` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(n_nodes) != 1L || !is.numeric(n_nodes) || n_nodes < 16) {
    stop("copsfm(): `n_nodes` must be a single number >= 16.", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(stats::terms(mf), mf)
  n <- length(y)
  k <- ncol(X)
  if (n <= k + 3L) {
    stop("copsfm(): ", n, " observations cannot identify ", k + 3L,
      " parameters.",
      call. = FALSE
    )
  }
  S <- if (isTRUE(inefdec)) 1 else -1
  cs <- .cop_spec(copula)

  ## Gauss-Legendre nodes on (0, 1), mapped to u in (0, Inf) by u = t/(1-t).
  ## The map is chosen so the integrand's mass -- which sits at moderate u --
  ## lands in the middle of the node range rather than in a tail.
  ##
  ## 128 nodes, not 64. Measured against the closed-form normal/half-normal
  ## density at independence, the maximum absolute error in log f is 1.3e-2 at
  ## 32 nodes, 2.5e-5 at 64 and 8.2e-14 at 128. 64 looks adequate and is not:
  ## an error of 1e-5 per observation is an error of 0.02 in a log-likelihood
  ## over 2000 points, which is the scale at which modes are being compared.
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  tt <- gl$nodes
  wt <- gl$weights
  uu <- tt / (1 - tt)
  jac <- 1 / (1 - tt)^2

  ## log f_eps(eps_i), the composed density, by quadrature over u.
  .log_dens <- function(eps, su, sv, cpar) {
    ## Rows are observations, columns quadrature nodes.
    U <- matrix(uu, nrow = length(eps), ncol = length(uu), byrow = TRUE)
    V <- eps + S * U
    lfv <- stats::dnorm(V, 0, sv, log = TRUE)
    lfu <- log(2) + stats::dnorm(U, 0, su, log = TRUE)
    lg <- lfv + lfu
    if (!identical(copula, "independent") && !is.null(cpar)) {
      ## Marginal CDFs, clamped off the endpoints: qnorm(0) is -Inf.
      w1 <- pmin(pmax(stats::pnorm(V / sv), 1e-12), 1 - 1e-12)
      w2 <- pmin(pmax(2 * stats::pnorm(U / su) - 1, 1e-12), 1 - 1e-12)
      lc <- .cop_logc(as.numeric(w1), as.numeric(w2), cpar, copula)
      if (any(!is.finite(lc))) return(rep(NA_real_, length(eps)))
      lg <- lg + matrix(lc, nrow = nrow(lg))
    }
    lg <- sweep(lg, 2, log(wt) + log(jac), "+")
    .log_row_sum_exp(lg)
  }

  ## Starting values: OLS plus Olson's moments, independence for the copula.
  ols <- stats::lm.fit(x = X, y = y)
  e0 <- as.numeric(ols$residuals)
  m2 <- mean(e0^2); m3 <- mean(e0^3)
  k3 <- sqrt(2 / pi) * (1 - 4 / pi)
  su0 <- if (S * m3 < 0) (S * m3 / k3)^(1 / 3) else 0.5 * stats::sd(e0)
  su0 <- if (is.finite(su0)) max(su0, 1e-3) else max(0.5 * stats::sd(e0), 1e-3)
  sv0 <- sqrt(max(m2 - (1 - 2 / pi) * su0^2, 1e-6))
  b0 <- as.numeric(ols$coefficients)
  if (attr(stats::terms(mf), "intercept") == 1L) {
    b0[1L] <- b0[1L] + S * su0 * sqrt(2 / pi)
  }
  start_v <- c(b0, log(su0), log(sv0), cs$par0)
  par_names <- c(colnames(X), "sigma_u", "sigma_v", cs$name)
  if (isTRUE(start_val)) names(start_v) <- par_names

  i_b <- seq_len(k); i_su <- k + 1L; i_sv <- k + 2L; i_c <- k + 3L

  like.fn <- function(th) {
    if (!all(is.finite(th))) return(cz$MAX_VALUE)
    su <- exp(pmin(th[i_su], 12)); sv <- exp(pmin(th[i_sv], 12))
    cpar <- th[i_c]
    if (cpar < cs$lo || cpar > cs$hi) return(cz$MAX_VALUE)
    eps <- S * as.numeric(y - X %*% th[i_b])
    ll <- .log_dens(eps, su, sv, cpar)
    if (any(!is.finite(ll))) return(cz$MAX_VALUE)
    ## The scaffold MINIMIZES; every likelihood here returns the negative sum.
    -sum(ll)
  }

  span <- pmax(10 * abs(b0), 10)
  lower1 <- c(b0 - span, log(1e-6), log(1e-6), cs$lo)
  upper1 <- c(b0 + span, log(1e4), log(1e4), cs$hi)

  Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v,
    lower.bobyqa = lower1, upper.bobyqa = upper1,
    maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, rhobeg = NA, rhoend = NA,
    verbose = verbose
  )
  start_v <- Opt.Bobyqa$start_v; bob1 <- Opt.Bobyqa$bob1

  Opt.Psoptim <- opt.psoptim(fn = like.fn, start_v,
    lower.psoptim = lower1, rand.psoptim = rand.psoptim,
    upper.psoptim = upper1, maxit.psoptim, psopt.TF = PSopt,
    rand.order = FALSE, verbose = verbose
  )
  start_v <- Opt.Psoptim$start_v; opt00 <- Opt.Psoptim$opt00

  Opt.Optim <- opt.optim(fn = like.fn, start_v = start_v,
    lower.optim = lower1, upper.optim = upper1, maxit.optim = maxit.optim,
    opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
  )
  start_v <- Opt.Optim$start_v; opt <- Opt.Optim$opt
  End.Time <- end.time(Start.Time)

  if (optHessian == FALSE & PSopt == FALSE) { opt <- bob1; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == FALSE & PSopt == TRUE)  { opt <- opt00; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == TRUE) {
    st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(diag(solve(opt$hessian))))
    }
  }

  th <- opt$par
  su <- exp(th[i_su]); sv <- exp(th[i_sv])
  par <- c(th[i_b], su, sv, th[i_c])
  se <- c(st_err[i_b], su * st_err[i_su], sv * st_err[i_sv], st_err[i_c])

  out <- matrix(NA_real_, 3L, length(par))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- par_names
  out[1, ] <- par; out[2, ] <- se; out[3, ] <- par / se

  ## E[u | eps] by the same quadrature the likelihood used, so the predictor and
  ## the objective cannot disagree about the model.
  eps <- S * as.numeric(y - X %*% th[i_b])
  Um <- matrix(uu, nrow = n, ncol = length(uu), byrow = TRUE)
  V <- eps + S * Um
  lg <- stats::dnorm(V, 0, sv, log = TRUE) + log(2) + stats::dnorm(Um, 0, su, log = TRUE)
  w1 <- pmin(pmax(stats::pnorm(V / sv), 1e-12), 1 - 1e-12)
  w2 <- pmin(pmax(2 * stats::pnorm(Um / su) - 1, 1e-12), 1 - 1e-12)
  lg <- lg + matrix(.cop_logc(as.numeric(w1), as.numeric(w2), th[i_c], copula), nrow = n)
  lg <- sweep(lg, 2, log(wt) + log(jac), "+")
  wts <- exp(lg - .log_row_sum_exp(lg))
  jlms <- as.numeric(rowSums(wts * Um))

  results <- list(
    t(out), c(opt), End.Time, start_v, "COP", formula, copula, th[i_c],
    jlms, exp(-jlms), S, n, as.integer(n_nodes),
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c(
    "out", "opt", "total_time", "start_v", "model_name", "formula", "copula",
    "copula_par", "jlms", "efficiency", "S", "nobs", "n_nodes",
    "coefficients", "std.errors", "t.values", "call"
  )
  results
}
