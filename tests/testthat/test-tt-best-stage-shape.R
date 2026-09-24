## Gap A50. ttsfm()'s TTHN branch compares the three optimizer stages and
## reports whichever attained the best objective. bobyqa() and psoptim() return
## their own shapes -- neither carries a hessian, and bobyqa() names its
## objective `fval` -- so substituting one raw left the branch calling
## colMeans(NULL) ("'x' must be an array of at least two dimensions") and left
## logLik() with no $value to read.

## optim()'s three stage-mates, in their real shapes.
.bob_like <- function(par, fval) {
  list(par = par, fval = fval, feval = 20L, ierr = 0L, msg = "normal")
}
.pso_like <- function(par, value) {
  list(par = par, value = value, counts = c(`function` = 20, iteration = 5),
       convergence = 0L, message = NULL)
}
.optim_like <- function(par, value) {
  list(par = par, value = value, counts = c(`function` = 20, gradient = 5),
       convergence = 0L, message = NULL, hessian = diag(2) * 3)
}

## The exact expression ttsfm.R runs on the chosen stage's hessian.
.se_expr <- function(h, np) {
  if (isTRUE(as.numeric(sum(colMeans(h))) == 0)) {
    rep(NA, np)
  } else {
    suppressWarnings(sqrt(diag(solve(h))))
  }
}

fq <- function(p) sum((p - c(1, 2))^2)

test_that("a winning bobyqa stage comes back in optim() shape", {
  cands <- list(bobyqa = .bob_like(c(1, 2), 0),
                psoptim = .pso_like(c(3, 3), 5),
                optim = .optim_like(c(4, 4), 13))
  bs <- .tt_best_stage(fq, cands, "TTHN")
  expect_identical(bs$which, "bobyqa")
  expect_true(is.matrix(bs$opt$hessian))
  expect_identical(dim(bs$opt$hessian), c(2L, 2L))
  ## bobyqa has no $value at all; before A50 logLik() read NULL and gave NA.
  expect_true(is.numeric(bs$opt$value))
  expect_equal(bs$opt$value, fq(c(1, 2)))
  ## The line that used to fail.
  expect_error(colMeans(bs$opt$hessian), NA)
  expect_error(.se_expr(bs$opt$hessian, 2), NA)
})

test_that("a winning psoptim stage comes back with a hessian", {
  cands <- list(bobyqa = .bob_like(c(4, 4), 13),
                psoptim = .pso_like(c(1, 2), 0),
                optim = .optim_like(c(3, 3), 5))
  bs <- .tt_best_stage(fq, cands, "TTHN")
  expect_identical(bs$which, "psoptim")
  expect_true(is.matrix(bs$opt$hessian))
  expect_error(.se_expr(bs$opt$hessian, 2), NA)
  ## Quadratic: the true hessian is 2I, so the SEs are computable and finite.
  expect_true(all(is.finite(.se_expr(bs$opt$hessian, 2))))
})

test_that("a winning optim stage is returned untouched", {
  op <- .optim_like(c(1, 2), 0)
  cands <- list(bobyqa = .bob_like(c(4, 4), 13),
                psoptim = .pso_like(c(3, 3), 5),
                optim = op)
  bs <- .tt_best_stage(fq, cands, "TTHN")
  expect_identical(bs$which, "optim")
  expect_identical(bs$opt, op)
})

test_that("an unusable numerical hessian degrades to NA, not to an error", {
  ## Defined at the candidate points and nowhere else, so numDeriv's
  ## perturbations throw and no finite hessian can be built.
  fragile <- function(p) {
    p <- as.numeric(p)
    if (isTRUE(all.equal(p, c(0.5, 0.5)))) return(1)
    if (isTRUE(all.equal(p, c(2, 2)))) return(5)
    stop("objective undefined here")
  }
  cands <- list(bobyqa = .bob_like(c(0.5, 0.5), 1),
                optim = .optim_like(c(2, 2), 5))
  bs <- .tt_best_stage(fragile, cands, "TTHN")
  expect_identical(bs$which, "bobyqa")
  expect_true(is.matrix(bs$opt$hessian))
  expect_true(all(is.na(bs$opt$hessian)))
  ## An all-NA hessian must reach NA standard errors without throwing.
  se <- .se_expr(bs$opt$hessian, 2)
  expect_error(.se_expr(bs$opt$hessian, 2), NA)
  expect_true(all(is.na(se)))
})

test_that("optHessian = FALSE skips the numerical hessian entirely", {
  cands <- list(bobyqa = .bob_like(c(1, 2), 0),
                optim = .optim_like(c(4, 4), 13))
  bs <- .tt_best_stage(fq, cands, "TTHN", optHessian = FALSE)
  expect_identical(bs$which, "bobyqa")
  expect_null(bs$opt$hessian)
  expect_equal(bs$opt$value, 0)
})

test_that("the substituted stage does not trigger the optim convergence warning", {
  cands <- list(bobyqa = .bob_like(c(1, 2), 0),
                optim = .optim_like(c(4, 4), 13))
  bs <- .tt_best_stage(fq, cands, "TTHN")
  ## ttsfm.R warns when opt$convergence is non-NULL and non-zero; a reported
  ## non-optim stage must not be described as a failed optim() run.
  expect_identical(bs$opt$convergence, 0L)
  expect_match(bs$opt$message, "not the best of the three stages")
})
