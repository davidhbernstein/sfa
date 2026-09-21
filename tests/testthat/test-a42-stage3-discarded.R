## A42. opt.optim() computes start_v/start_feval as the better of the point it
## was handed and optim()'s own result, but every call site reports `opt`, so
## whenever optim() ended above its start the better answer was discarded.
## psfm(model_name = "TRE_Z") did this on 14 of 50 simulated fits, by up to
## 1.12 log-likelihood points.
##
## Reporting start_v unconditionally is the fix that was tried and reverted:
## in the NR spike case start_feval > opt$value is FALSE, so start_v still
## holds the spike. test-ttsfm-numerics.R pins that. Hence the credibility
## gate, and hence the two tests below that demand opposite answers from it.

test_that(".sfa_point_credible rejects a collapse and accepts an interior point", {
  ## Contributions must VARY: a constant per-observation vector is itself one
  ## of the degeneracies the gate rejects, so a flat fixture would test the
  ## wrong branch.
  fn <- function(p, per_obs = FALSE) {
    v <- -(p[1]^2 + seq_len(5))
    if (per_obs) v else -sum(v)
  }
  lo <- c(0, 0)
  up <- c(10, 10)

  ## Ordinary interior point.
  expect_true(.sfa_point_credible(fn, c(2, 3), 5, lo, up))

  ## A summed log-likelihood of exactly zero means every density evaluated to
  ## exactly 1 -- the A49 collapse. Negative zero must be caught too.
  expect_false(.sfa_point_credible(fn, c(2, 3), 0, lo, up))
  expect_false(.sfa_point_credible(fn, c(2, 3), -0, lo, up))

  expect_false(.sfa_point_credible(fn, c(2, 3), NA_real_, lo, up))
  expect_false(.sfa_point_credible(fn, c(2, 3), Inf, lo, up))

  ## On a bound: the signature of a scale driven to its floor, which is what
  ## the spike the guard must not take looks like.
  expect_false(.sfa_point_credible(fn, c(0, 3), 5, lo, up))
  expect_false(.sfa_point_credible(fn, c(2, 10), 5, lo, up))

  ## Every per-observation contribution identical, or all exactly zero.
  fn_zero <- function(p, per_obs = FALSE) {
    v <- rep(0, 5)
    if (per_obs) v else -sum(v)
  }
  expect_false(.sfa_point_credible(fn_zero, c(2, 3), 1, lo, up))
})

test_that("a closure with no per_obs argument is still usable", {
  ## Most of the scaffold's callers pass a plain objective. The gate must fall
  ## back to the cheap checks rather than erroring on the missing argument.
  fn <- function(p) sum(p^2)
  expect_true(.sfa_point_credible(fn, c(1, 1), 2, c(-5, -5), c(5, 5)))
})

test_that("opt.optim keeps a better stage-2 point instead of discarding it", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(
    t = 5, N = 40, rand = 2008, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
    sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5
  ))
  fit <- suppressWarnings(psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
    data = d, individual = "name"
  ))
  ## Before the guard this reported 294.1174926; stage 2 had reached
  ## 293.9555300 and optim() ended above it.
  expect_lt(as.numeric(fit$opt$value), 294.0)
  expect_equal(as.numeric(fit$opt$value), 293.95553, tolerance = 1e-4)

  ## What is reported must stay self-consistent: logLik() reads $opt$value, and
  ## the Hessian must describe the point being reported, not the one optim
  ## wandered to.
  expect_equal(-as.numeric(stats::logLik(fit)), as.numeric(fit$opt$value),
    tolerance = 1e-9
  )
  expect_true(all(is.finite(fit$opt$hessian)))
  expect_true(any(is.finite(as.numeric(fit$std.errors))))
})

test_that("the NR spike is still escaped rather than reported", {
  skip_on_cran()
  ## The counterexample the guard must not break -- same data as
  ## test-ttsfm-numerics.R. If the gate ever accepts the spike, logLik jumps to
  ## +3.6e17 instead of the hundreds.
  set.seed(3)
  n <- 400
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  d <- data.frame(
    y = 1 + 0.5 * x1 + 0.5 * x2 + stats::rnorm(n, 0, 0.4) - stats::rexp(n, 1),
    x1 = x1, x2 = x2
  )
  f <- suppressWarnings(sfm(y ~ x1 + x2, model_name = "NR", data = d))
  ll <- as.numeric(stats::logLik(f))
  expect_true(is.finite(ll))
  expect_lt(ll, 0)
  expect_gt(ll, -1e4)
})
