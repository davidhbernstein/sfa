## A53. A stage substituted for optim()'s result must report a Hessian from the
## SAME estimator optim(hessian = TRUE) uses, and at a tie must carry optim()'s
## own Hessian across unchanged. Otherwise the standard errors depend on which
## stage won a tie. These tests pin the estimator, not where any fit lands.

.a53_nll <- local({
  set.seed(53)
  n <- 150
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n, 0, 0.4)
  function(p) -sum(stats::dnorm(y, p[1] + p[2] * x, exp(p[3]), log = TRUE))
})
.a53_par <- c(0.9, 0.55, log(0.45))

test_that("optimHess() is the estimator optim(hessian = TRUE) reports", {
  ## The premise of the fix: at optim()'s own point, optimHess() with default
  ## control reproduces optim()'s stored Hessian exactly. opt.optim() passes no
  ## parscale or ndeps, so the defaults are the ones that matter.
  op <- stats::optim(.a53_par, .a53_nll, method = "L-BFGS-B",
                     lower = c(-5, -5, -5), upper = c(5, 5, 5), hessian = TRUE)
  expect_identical(stats::optimHess(op$par, .a53_nll), op$hessian)
})

test_that(".sfa_substitute_hessian away from a tie equals optimHess()", {
  h <- .sfa_substitute_hessian(.a53_nll, .a53_par)
  expect_identical(h, stats::optimHess(.a53_par, .a53_nll))

  ## A reference point that is NOT a tie is ignored.
  ref <- list(par = .a53_par + 1e-4, hessian = diag(3))
  expect_identical(.sfa_substitute_hessian(.a53_nll, .a53_par, ref = ref),
                   stats::optimHess(.a53_par, .a53_nll))
})

test_that("at a tie the reference Hessian is carried across as-is", {
  carried <- matrix(c(4, 1, 0, 1, 3, 0, 0, 0, 2), 3, 3)
  ## 1e-10 relative: the size of the recorded TTHN tie (8.87e-11).
  ref <- list(par = .a53_par * (1 + 1e-10), hessian = carried)
  expect_identical(.sfa_substitute_hessian(.a53_nll, .a53_par, ref = ref), carried)

  ## Just outside the tolerance it is recomputed instead.
  ref_far <- list(par = .a53_par + 2 * .SFA_HESS_TIE_TOL * pmax(1, abs(.a53_par)),
                  hessian = carried)
  expect_identical(.sfa_substitute_hessian(.a53_nll, .a53_par, ref = ref_far),
                   stats::optimHess(.a53_par, .a53_nll))

  ## A non-finite reference Hessian is not carried.
  ref_na <- list(par = .a53_par, hessian = matrix(NA_real_, 3, 3))
  expect_identical(.sfa_substitute_hessian(.a53_nll, .a53_par, ref = ref_na),
                   stats::optimHess(.a53_par, .a53_nll))
})

test_that(".tt_as_optim uses optimHess(), not numDeriv", {
  z <- list(par = .a53_par, fval = .a53_nll(.a53_par))
  out <- .tt_as_optim(z, .a53_nll, .a53_nll(.a53_par), "bobyqa", TRUE)
  expect_identical(out$hessian, stats::optimHess(.a53_par, .a53_nll))
})

test_that(".tt_best_stage carries optim()'s Hessian across a tie", {
  carried <- stats::optimHess(.a53_par, .a53_nll) + 0.25
  ## bobyqa wins by more than optim()'s reltol margin, at a point 1e-10 away.
  bump <- function(p) .a53_nll(p) + if (isTRUE(all.equal(p, .a53_par, tolerance = 0))) 0 else 1e-3
  cands <- list(
    bobyqa = list(par = .a53_par, fval = 0, feval = 1L, ierr = 0L, msg = ""),
    optim = list(par = .a53_par * (1 + 1e-10), value = 0, counts = c(1L, 1L),
                 convergence = 0L, message = NULL, hessian = carried)
  )
  bs <- .tt_best_stage(bump, cands, "TTHN")
  expect_identical(bs$which, "bobyqa")
  expect_identical(bs$opt$par, .a53_par)
  expect_identical(bs$opt$hessian, carried)
})

test_that(".tt_best_stage away from a tie computes optimHess() at the winner", {
  far <- .a53_par + c(0.5, 0.5, 0.5)
  cands <- list(
    bobyqa = list(par = .a53_par, fval = 0, feval = 1L, ierr = 0L, msg = ""),
    optim = list(par = far, value = 0, counts = c(1L, 1L), convergence = 0L,
                 message = NULL, hessian = diag(3))
  )
  bs <- .tt_best_stage(.a53_nll, cands, "TTHN")
  expect_identical(bs$which, "bobyqa")
  expect_identical(bs$opt$hessian, stats::optimHess(.a53_par, .a53_nll))
})
