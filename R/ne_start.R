## Starting values for the normal-exponential frontier (Bernstein, Parmeter and Wright,
## "Starting Values for the Normal-Exponential Stochastic Frontier Model"): four estimators
## of sigma_u from OLS residuals; the default `bc` is mu1- corrected by h(lambda)^w.
## See notes/code_history/ne_start.md.

## h(lambda) = plim(mu1-)/sigma_u >= 1, the exact asymptotic bias factor of the mu1-
## start, formed in logs: its exp overflows below lambda ~ 0.03 as its Phi underflows.
.ne_bias_factor <- function(lambda) {
  lambda <- pmax(lambda, 1e-8)
  Psi <- exp(1 / (2 * lambda^2) - 1 + stats::pnorm(lambda - 1 / lambda, log.p = TRUE))
  (Psi + stats::dnorm(lambda) / lambda) / (stats::pnorm(-lambda) + Psi)
}

## N * asymptotic variance of each start, in units of sigma_u^2.  Used by the
## `adaptive` rule and by the theory tables; not needed by `bc`.
.ne_avar_mu1 <- function(lambda) {
  lambda <- pmax(lambda, 1e-8)
  Psi <- exp(1 / (2 * lambda^2) - 1 + stats::pnorm(lambda - 1 / lambda, log.p = TRUE))
  p <- stats::pnorm(-lambda) + Psi
  A <- Psi + stats::dnorm(lambda) / lambda
  M2 <- 2 * Psi + (1 + lambda^-2) * stats::pnorm(-lambda) + stats::dnorm(lambda) / lambda
  Vm <- M2 - A^2 / p
  ## g: first-order effect of the OLS intercept shift on mu1-.  Omitting it
  ## overstates the variance by ~25% at lambda = 0.5.
  g <- 1 - (A / p) * Psi / p
  (1 / p^2) * Vm + g^2 * (1 + lambda^-2) - (2 * g / p) * Vm
}

.ne_avar_mu3 <- function(lambda) {
  s <- 1 + pmax(lambda, 1e-8)^-2
  (26 + 9 * s + s^3) / 6
}

## lambda*(N), where the mu1- and mu3 starts have equal first-order MSE. Memoised on N:
## a Monte Carlo asks for the same few sample sizes thousands of times.
.ne_lambda_star_cache <- new.env(parent = emptyenv())

.ne_lambda_star <- function(N, interval = c(0.3, 40)) {
  vapply(N, function(n) {
    key <- paste0("n", n, "_", interval[1], "_", interval[2])
    if (!is.null(hit <- .ne_lambda_star_cache[[key]])) {
      return(hit)
    }
    val <- .ne_lambda_star_one(n, interval)
    assign(key, val, envir = .ne_lambda_star_cache)
    val
  }, numeric(1))
}

.ne_lambda_star_one <- function(N, interval) {
  vapply(N, function(n) {
    ff <- function(l) (.ne_bias_factor(l) - 1)^2 + (.ne_avar_mu1(l) - .ne_avar_mu3(l)) / n
    if (!is.finite(ff(interval[1])) || !is.finite(ff(interval[2])) ||
      ff(interval[1]) * ff(interval[2]) > 0) {
      return(NA_real_)
    }
    stats::uniroot(ff, interval = interval, tol = 1e-10)$root
  }, numeric(1))
}

## All four sigma_u starts from a vector of OLS residuals, plus the matching
## sigma_v and the intercept correction E[u] = sigma_u.
##
## `rule` selects which one the caller gets back in $sigma_u.  Returns every
## candidate in $candidates either way, so a caller can compare them.
.ne_moment_start <- function(epsilon_hat, rule = c("bc", "min", "mu1", "mu3", "adaptive"),
                             shrink = TRUE) {
  rule <- match.arg(rule)
  e <- as.numeric(epsilon_hat)
  e <- e[is.finite(e)]
  n <- length(e)
  eps_floor <- .Machine$double.eps

  if (n < 4L) {
    return(list(sigma_u = 0.1, sigma_v = 0.1, eu = 0.1, lambda = 1,
                rule = rule, wrong_skew = NA, candidates = NULL))
  }

  e <- e - mean(e) ## OLS residuals are already centred; harmless if not
  m2 <- mean(e^2)
  m3 <- mean(e^3)

  ## mu1-: always defined and always strictly positive, which is exactly why
  ## it is the fallback when the third moment has the wrong sign.
  neg <- e < 0
  su_1 <- if (any(neg)) max(-mean(e[neg]), eps_floor) else eps_floor

  ## mu3: the COLS inversion.  Undefined when the residuals are skewed the
  ## wrong way, which happens with non-trivial probability at low lambda --
  ## about 17% at lambda = 0.75, N = 100.
  wrong_skew <- !is.finite(m3) || m3 >= 0
  su_3 <- if (wrong_skew) NA_real_ else max((-m3 / 2)^(1 / 3), eps_floor)

  su_min <- if (wrong_skew) su_1 else min(su_1, su_3)

  ## lambda-hat comes from the mu3 start because that one is consistent.
  ## Under wrong skew there is no usable lambda-hat; treat it as high-noise,
  ## which is the regime that produced the wrong skew in the first place.
  lam_hat <- if (wrong_skew) NA_real_ else su_3 / sqrt(max(m2 - su_3^2, eps_floor))

  ## bc: divide out the exact asymptotic bias, shrunk by the share of MSE it
  ## explains.  With no lambda-hat, fall back to mu1- uncorrected rather than
  ## to an arbitrary correction.
  su_bc <- if (is.na(lam_hat)) {
    su_1
  } else {
    hh <- .ne_bias_factor(lam_hat)
    w <- if (shrink) {
      bb2 <- (hh - 1)^2
      vv <- .ne_avar_mu1(lam_hat) / n
      if (bb2 + vv <= 0) 0 else bb2 / (bb2 + vv)
    } else {
      1
    }
    max(su_1 / hh^w, eps_floor)
  }

  ## adaptive: pick the start with the smaller predicted MSE at lambda-hat.
  ls <- .ne_lambda_star(n)
  su_ad <- if (is.na(lam_hat) || is.na(ls)) su_1 else if (lam_hat > ls) su_1 else su_3

  cand <- c(mu1 = su_1, mu3 = su_3, min = su_min, bc = su_bc, adaptive = su_ad)
  su <- unname(cand[rule])
  if (!is.finite(su) || su <= 0) su <- su_1

  ## sigma_v from the second moment: Var(eps) = sigma_u^2 + sigma_v^2.
  sv <- sqrt(max(m2 - su^2, eps_floor))

  list(sigma_u = su, sigma_v = sv, eu = su, lambda = su / sv,
       rule = rule, wrong_skew = wrong_skew, lambda_hat = lam_hat,
       lambda_star = ls, candidates = cand, moments = c(m2 = m2, m3 = m3))
}

## Drop-in replacement for the `model_name == "NE"` branch of start_cs().
## Returns start_v in the package's NE order: (sigma_v, sigma_u, beta_0,
## beta_slopes).  Only the intercept is corrected -- the OLS slopes are
## already consistent.
.ne_start <- function(epsilon_hat, beta_0_st, beta_hat, rule = "bc", shrink = TRUE) {
  ms <- .ne_moment_start(epsilon_hat, rule = rule, shrink = shrink)
  if (is.na(beta_0_st)) {
    unname(c(ms$sigma_v, ms$sigma_u, beta_hat))
  } else {
    unname(c(ms$sigma_v, ms$sigma_u, unname(beta_0_st) + ms$eu, beta_hat))
  }
}
