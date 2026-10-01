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
  ## rand 2, 4, 5, 11 and 12 reach the boundary on this DGP; 1, 3, 6-10 do not.
  boundary <- c(2L, 4L, 5L, 11L, 12L)
  interior <- c(1L, 3L, 6L, 7L, 8L, 9L, 10L)

  ratio_of <- function(f) {
    p <- f$out[, "par"]
    lam <- p[["lambda"]]; sg <- p[["sigma"]]
    su <- lam * sg / sqrt(1 + lam^2); sv <- su / lam
    sv / sqrt(su^2 + sv^2)
  }
  fit <- function(r) {
    d <- data_gen_cs(N = 150, rand = r, cons = .5, beta1 = .5, beta2 = .5,
      sig_u = 1, sig_v = .5, mu = .5, a = 5)
    warned <- FALSE
    f <- withCallingHandlers(
      sfm(y_pcs_tn ~ x1 + x2, data = d, model_name = "NTN"),
      warning = function(w) {
        if (grepl("sigma_v has collapsed", conditionMessage(w))) warned <<- TRUE
        invokeRestart("muffleWarning")
      })
    list(f = f, warned = warned)
  }

  for (r in boundary) {
    z <- fit(r)
    expect_true(isTRUE(z$f$ntn_sigma_v_at_bound), info = paste("rand", r))
    expect_true(z$warned, info = paste("rand", r))
    ## The flag has to mean what it says: no standard errors survive here.
    expect_true(all(is.na(z$f$out[, "st_err"])), info = paste("rand", r))
  }
  for (r in interior) {
    z <- fit(r)
    expect_false(isTRUE(z$f$ntn_sigma_v_at_bound), info = paste("rand", r))
    expect_false(z$warned, info = paste("rand", r))
    ## ... and where it is not set, the fit is ordinary.
    expect_true(all(is.finite(z$f$out[, "st_err"])), info = paste("rand", r))
  }
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
