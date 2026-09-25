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

test_that("the substituted hessian is on the same scale optim() would have used", {
  ## The point of the substitution is that SEs stay comparable no matter which
  ## stage won. If numDeriv's hessian disagreed with optim(hessian = TRUE)'s,
  ## a reported SE would mean different things on different seeds.
  set.seed(11)
  n <- 200
  x <- rnorm(n)
  y <- 1.4 + 0.7 * x + rnorm(n, 0, 0.5)
  ## A real negative log-likelihood, the shape every ttsfm closure returns.
  nll <- function(p) {
    -sum(stats::dnorm(y, p[1] + p[2] * x, exp(p[3]), log = TRUE))
  }
  op <- stats::optim(c(1, 1, -1), nll, method = "BFGS", hessian = TRUE)

  cands <- list(bobyqa = .bob_like(op$par, op$value),
                optim = .optim_like(c(9, 9, 9), nll(c(9, 9, 9))))
  bs <- .tt_best_stage(nll, cands, "TTHN")
  expect_identical(bs$which, "bobyqa")

  se_sub <- suppressWarnings(sqrt(diag(solve(bs$opt$hessian))))
  se_opt <- suppressWarnings(sqrt(diag(solve(op$hessian))))
  expect_equal(se_sub, se_opt, tolerance = 1e-5)
  ## And they are real standard errors, not an artefact: the slope SE of an
  ## OLS fit with the same design agrees.
  expect_equal(se_sub[2], unname(summary(stats::lm(y ~ x))$coefficients[2, 2]),
               tolerance = 5e-3)
})
