## NTN's sigma_v boundary, which used to be reported by silence.
##
## On some samples NTN has no interior maximum in lambda: the profile
## log-likelihood rises monotonically toward a finite limit as lambda -> Inf.
## Measured by continuation from the reported solution on rand = 2 --
## -190.3374 at lambda = 1e1, -188.5223 at 1e2, -188.1011 at 1e5, -188.1003 at
## the fitted 9.72e5, and flat at -188.1001876 from 1e9 through 1e15. So the
## boundary IS the supremum and running to lambda = 1e6 is correct; but the
## likelihood is flat in lambda there, the Hessian is singular, and every
## standard error comes back NA.
##
## Before this, the user got lambda = 972468, six NA standard errors, and no
## message of any kind. tHN's mirror case (sigma_u -> 0) has warned and set
## $thn_sigma_u_at_bound since 1.2.0; this is the same treatment for the other
## component.

test_that("NTN flags and warns when sigma_v collapses, and only then", {
  skip_on_cran()
  ## Which seeds reach the corner is where the optimizer lands, and it moved
  ## when the ladder learned to escape the mu -> -Inf ridge (rand 9 and 10).
  ## So instead of a list of seeds, check each fit against the exact
  ## lambda -> Inf limit: a deterministic frontier with u = -eps truncated
  ## normal, maximised directly.
  exact_limit <- function(d, start) {
    X <- cbind(1, d$x1, d$x2)
    y <- d$y_pcs_tn
    ## th = (slopes, mu, log sigma_u, log offset): the intercept sits
    ## exp(offset) above the frontier, so every u = -eps stays >= 0.
    h <- function(th) {
      b <- c(0, th[1:2])
      b[1] <- max(y - X %*% b) + exp(th[5])
      u <- -(y - X %*% b)
      s <- exp(th[4])
      -sum(dnorm((u - th[3]) / s, log = TRUE) - log(s) -
        pnorm(th[3] / s, log.p = TRUE))
    }
    best <- -Inf
    for (st in list(c(start, -10), c(.5, .5, .5, 0, -10))) {
      o <- stats::nlminb(st, h)
      o <- stats::optim(o$par, h, method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-12))
      best <- max(best, -o$value)
    }
    best
  }

  n_flag <- 0L
  for (r in 1:12) {
    d <- data_gen_cs(N = 150, rand = r, cons = .5, beta1 = .5, beta2 = .5,
      sig_u = 1, sig_v = .5, mu = .5, a = 5)
    warned <- FALSE
    f <- withCallingHandlers(
      sfm(y_pcs_tn ~ x1 + x2, data = d, model_name = "NTN"),
      warning = function(w) {
        if (grepl("sigma_v has collapsed", conditionMessage(w))) warned <<- TRUE
        invokeRestart("muffleWarning")
      })
    p <- f$out[, "par"]
    su <- p[["lambda"]] * p[["sigma"]] / sqrt(1 + p[["lambda"]]^2)
    lim <- exact_limit(d, c(p[c("x1", "x2")], p[["mu"]], log(su)))
    ll <- as.numeric(logLik(f))
    flag <- isTRUE(f$ntn_sigma_v_at_bound)
    n_flag <- n_flag + flag
    info <- paste("rand", r)
    ## The flag, the warning and the missing SEs travel together.
    expect_identical(warned, flag, info = info)
    if (flag) {
      expect_true(all(is.na(f$out[, "st_err"])), info = info)
      ## At the corner the ladder's top rung (lambda ~ 1e11) ends in nlminb
      ## "false convergence"; Ubuntu oldrel stopped 0.0014 below the exact
      ## limit on rand 9, macOS within 1e-5.
      expect_lt(abs(ll - lim), 5e-3, label = info)
    } else {
      expect_true(all(is.finite(f$out[, "st_err"])), info = info)
      ## Not flagged means the corner is not higher than the fit.
      expect_gte(ll, lim - 1e-4, label = info)
    }
  }
  ## Both regimes occur on this DGP, so neither branch is vacuous.
  expect_gte(n_flag, 4L)
  expect_lte(n_flag, 10L)
})

test_that("the boundary criterion is not a tuned cut", {
  skip_on_cran()
  ## The threshold is sigma_v/sigma < 1e-3. What makes that safe is not the
  ## number but the gap: measured over these 12 samples the boundary fits sit
  ## at 4.9e-08 to 1.3e-06 and the interior fits at 0.126 to 0.442 -- five
  ## orders of magnitude apart, with the threshold in the middle of empty
  ## space. Asserted so that a future change which narrows that gap fails here
  ## rather than silently making the flag sensitive to the cut.
  lo <- numeric(0); hi <- numeric(0)
  for (r in 1:12) {
    d <- data_gen_cs(N = 150, rand = r, cons = .5, beta1 = .5, beta2 = .5,
      sig_u = 1, sig_v = .5, mu = .5, a = 5)
    f <- suppressWarnings(sfm(y_pcs_tn ~ x1 + x2, data = d, model_name = "NTN"))
    p <- f$out[, "par"]
    lam <- p[["lambda"]]; sg <- p[["sigma"]]
    su <- lam * sg / sqrt(1 + lam^2); sv <- su / lam
    rt <- sv / sqrt(su^2 + sv^2)
    if (isTRUE(f$ntn_sigma_v_at_bound)) lo <- c(lo, rt) else hi <- c(hi, rt)
  }
  expect_gt(length(lo), 0)
  expect_gt(length(hi), 0)
  ## At least two orders of clearance on each side of 1e-3.
  expect_lt(max(lo), 1e-5)
  expect_gt(min(hi), 1e-1)
})
