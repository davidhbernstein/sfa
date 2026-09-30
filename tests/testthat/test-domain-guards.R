## Out-of-domain penalties in the likelihood closures.
##
## Every branch guards against a non-positive scale, but the penalty has to be
## a large FINITE number. optim() differences the objective to get a gradient,
## and differencing .Machine$double.xmax overflows to a non-finite value, which
## aborts the fit with "non-finite finite-difference value" rather than steering
## the search away. NGE used xmax and died on 3 of 45 fits at N = 150 -- including
## at sigma_u = 1, sigma_v = 0.3 -- before this was changed. NU and NE carried
## the same construct.

test_that("no likelihood branch returns .Machine$double.xmax as its penalty", {
  ## Source-level guard, and only meaningful when the sources are on disk --
  ## under R CMD check the package is installed and R/ is gone, so this skips.
  ## readLines() on a missing path ERRORS rather than returning empty, so the
  ## existence check has to come first.
  ## Every file carrying a likelihood closure, not just sfm.R. Checking only
  ## sfm.R is how ttsfm.R kept SEVEN of these until 2026-09-04: TTNE, TTHN and
  ## TTNLS each returned xmax, and TTHN's fits were failing with optim()'s
  ## ABNORMAL_TERMINATION_IN_LNSRCH -- the same non-finite-gradient mechanism
  ## this test exists to prevent, in a file it never read.
  want <- c("sfm.R", "ttsfm.R", "zsfm.R", "lcsfm.R", "psfm.R", "selsfm.R",
    "ivsfm.R", "copsfm.R")
  roots <- c("../../R", "../../../R", "R")
  root <- roots[dir.exists(roots)]
  skip_if(!length(root), "R/ source not reachable (installed package)")
  files <- file.path(root[1], want)
  files <- files[file.exists(files)]
  skip_if(!length(files), "no likelihood sources found")
  src <- unlist(lapply(files, readLines, warn = FALSE))
  ## `return(.Machine$double.xmax)` is the exact construct that overflows
  ## optim()'s finite-difference gradient.
  bad <- grep("return\\(\\.Machine\\$double\\.xmax\\)", src, value = TRUE)
  expect_equal(length(bad), 0)
})

## WIDENED 2026-09-30 (gap A56a). The test above matches only the literal
## `return(.Machine$double.xmax)`, and sfm()'s tHN branch wrote the penalty as
##
##     rep(-.Machine$double.xmax / length(eps), length(eps))
##
## which is not a `return(...)` and so walked straight through the one test
## written to prevent it. Summed over n that IS double.xmax, so it is strictly
## worse than the construct the test did catch: the first finite difference
## gives Inf, not merely a huge number. It cost 25 of 84 tHN fits every
## standard error.
##
## The distinction below is mechanical, not stylistic:
##   bare      xmax / n  summed = xmax (1.8e308)      -> differences to Inf
##   sqrt(xmax / n)      summed = sqrt(xmax * n) ~ 2e155 -> huge but FINITE
## The second form is the general `!is.finite(like)` guard, which still lives in
## sfm.R, zsfm.R, lcsfm.R and ttsfm.R and is tracked as the still-open gap A56.
## It is deliberately NOT failed here; this test draws the line at the form that
## is non-finite on the first difference.
test_that("no likelihood uses .Machine$double.xmax as a penalty VALUE", {
  want <- c("sfm.R", "ttsfm.R", "zsfm.R", "lcsfm.R", "psfm.R", "selsfm.R",
    "ivsfm.R", "copsfm.R")
  roots <- c("../../R", "../../../R", "R")
  root <- roots[dir.exists(roots)]
  skip_if(!length(root), "R/ source not reachable (installed package)")
  files <- file.path(root[1], want)
  files <- files[file.exists(files)]
  skip_if(!length(files), "no likelihood sources found")

  bad <- character(0)
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    for (i in seq_along(ln)) {
      s <- sub("^[[:space:]]+", "", ln[i])
      if (startsWith(s, "#")) next                       # comment line
      if (!grepl("\\.Machine\\$double\\.xmax", s)) next
      ## Allowed: wrapped in sqrt(), i.e. the A56 general guard.
      if (grepl("sqrt\\([^()]*\\.Machine\\$double\\.xmax", s)) next
      ## Allowed: the constants table's own MAX_VALUE alias.
      if (grepl("MAX_VALUE", s)) next
      bad <- c(bad, sprintf("%s:%d: %s", basename(f), i, s))
    }
  }
  expect_equal(bad, character(0),
    info = paste0(
      "A bare .Machine$double.xmax used as a penalty value sums to xmax over ",
      "the observations and makes optim()'s finite difference non-finite. Use ",
      ".SFA_CONSTANTS$DOMAIN_PENALTY_PER_OBS. Offending lines:\n",
      paste(bad, collapse = "\n")))
})

test_that("the tHN domain penalty is the swept per-observation constant", {
  ## Pins the VALUE as well as the form: A56a measured standard-error recovery
  ## against the summed penalty (1e12 recovers nothing, 1e6 and below recover
  ## all five test fits), and per observation rather than as a fixed total so
  ## the margin over a legitimate optimum does not erode as n grows.
  expect_true(is.finite(.SFA_CONSTANTS$DOMAIN_PENALTY_PER_OBS))
  expect_gt(.SFA_CONSTANTS$DOMAIN_PENALTY_PER_OBS, 1e3)
  expect_lt(.SFA_CONSTANTS$DOMAIN_PENALTY_PER_OBS, 1e6)
  ## Summed at a realistic n it must stay well inside the band where the
  ## Hessian is still computable (measured: 1e8 total already costs one fit
  ## of five its standard errors).
  expect_lt(.SFA_CONSTANTS$DOMAIN_PENALTY_PER_OBS * 200, 1e8)
})

test_that("models with domain guards fit without aborting", {
  skip_on_cran()
  ## These configurations produced "non-finite finite-difference value" before
  ## the penalties were made finite.
  cfgs <- list(c(0.3, 0.4), c(1, 0.3), c(0.2, 0.8))
  for (mn in c("NGE", "NU", "NE")) {
    yv <- switch(mn, NGE = "y_pcs_ge", NU = "y_pcs_u", NE = "y_pcs_e")
    errs <- 0L
    for (cfg in cfgs) {
      for (r in 1:5) {
        d <- data_gen_cs(N = 150, rand = r, sig_u = cfg[1], sig_v = cfg[2], cons = 0.5,
                         beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
        f <- try(suppressWarnings(suppressMessages(
          sfm(as.formula(paste(yv, "~ x1 + x2")), model_name = mn, data = d))), silent = TRUE)
        if (inherits(f, "try-error")) errs <- errs + 1L
      }
    }
    expect_equal(errs, 0L, info = mn)
  }
})
