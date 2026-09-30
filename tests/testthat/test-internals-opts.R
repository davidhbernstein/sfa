## Direct tests of the optimizer-stage wrappers and small helpers that every
## entry point runs through but that had no test of their own. Each
## expectation is derived from the helper's contract -- a stage never hands
## back a worse point than it was given, a switched-off stage is a no-op, the
## caller's random stream survives -- not copied from current output.

fq <- function(p) sum((p - c(1, -2))^2) + 5      # minimum 5 at (1, -2)

test_that("opt.nlminb(): off is a no-op, on reaches the optimum it reports", {
  off <- opt.nlminb(fq, c(0, 0), lower.nlminb = -Inf, nlminb.TF = FALSE)
  expect_identical(off$start_v, c(0, 0))
  expect_identical(off$start_feval, fq(c(0, 0)))
  expect_null(off$nlm1)

  on <- opt.nlminb(fq, c(0, 0), lower.nlminb = -Inf, nlminb.TF = TRUE)
  expect_equal(on$start_v, c(1, -2), tolerance = 1e-6)
  ## The reported value is the objective AT the reported point.
  expect_equal(on$start_feval, fq(on$start_v))
  expect_identical(on$start_v, on$nlm1$par)
})

test_that("opt.nlminb(): an objective that fails keeps the start", {
  ## nlminb() stops on a non-finite value; the stage must cost the stage, not
  ## the fit.
  fbad <- function(p) if (identical(p, c(0, 0))) 3 else NA_real_
  r <- suppressWarnings(opt.nlminb(fbad, c(0, 0), lower.nlminb = -Inf, nlminb.TF = TRUE))
  expect_identical(r$start_v, c(0, 0))
  expect_identical(r$start_feval, 3)
})

test_that("opt.bobyqa(): off is a no-op; on never returns a worse point", {
  off <- opt.bobyqa(fq, c(0, 0), lower.bobyqa = c(-5, -5), upper.bobyqa = c(5, 5),
    maxit.bobyqa = 200, bob.TF = FALSE, verbose = FALSE)
  expect_identical(off$start_v, c(0, 0))
  expect_null(off$bob1)

  on <- opt.bobyqa(fq, c(0, 0), lower.bobyqa = c(-5, -5), upper.bobyqa = c(5, 5),
    maxit.bobyqa = 2000, bob.TF = TRUE, rhobeg = 0.5, rhoend = 1e-8,
    verbose = FALSE)
  expect_equal(on$start_v, c(1, -2), tolerance = 1e-5)
  expect_equal(on$start_feval, fq(on$start_v))
  expect_lte(on$start_feval, fq(c(0, 0)))
})

test_that("opt.bobyqa(): a start outside the bounds is moved just inside (A35)", {
  ## bobyqa() itself stops on such a start. Only the offending coordinate
  ## moves, and only by a relative 1e-6.
  r <- opt.bobyqa(fq, c(-1e-9, 0), lower.bobyqa = c(1e-7, -5), upper.bobyqa = c(5, 5),
    maxit.bobyqa = 2000, bob.TF = TRUE, rhobeg = 0.5, rhoend = 1e-8,
    verbose = FALSE)
  expect_true(all(r$start_v >= c(1e-7, -5)))
  expect_equal(r$start_v, c(1, -2), tolerance = 1e-5)
})

test_that("opt.psoptim(): off is a no-op; a seeded run is reproducible and leaves the caller's stream", {
  off <- opt.psoptim(fq, c(0, 0), lower.psoptim = c(-5, -5), upper.psoptim = c(5, 5),
    maxit.psoptim = 20, psopt.TF = FALSE, verbose = FALSE, rand.psoptim = 1)
  expect_identical(off$start_v, c(0, 0))
  expect_null(off$opt00)

  set.seed(42)
  s0 <- .Random.seed
  run <- function() opt.psoptim(fq, c(0, 0), lower.psoptim = c(-5, -5),
    upper.psoptim = c(5, 5), maxit.psoptim = 30, psopt.TF = TRUE,
    rand.order = FALSE, verbose = FALSE, rand.psoptim = 7)
  a <- run()
  expect_identical(.Random.seed, s0)
  b <- run()
  expect_identical(a$start_v, b$start_v)
  expect_lte(a$start_feval, fq(c(0, 0)))
  expect_equal(a$start_feval, fq(a$start_v))
})

test_that(".rng_snapshot()/.rng_restore() round-trip, including 'no stream'", {
  set.seed(3)
  s <- .rng_snapshot()
  runif(5)
  .rng_restore(s)
  expect_identical(.Random.seed, s)

  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  rm(".Random.seed", envir = globalenv())
  s0 <- .rng_snapshot()
  expect_null(s0)
  runif(1)
  .rng_restore(s0)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that(".safe_symmetrize() averages with the transpose and drops only spurious dimnames", {
  M <- matrix(c(2, 1, 3, 4), 2)
  S <- .safe_symmetrize(M)
  expect_identical(S, (M + t(M)) / 2)
  expect_true(isSymmetric(S))

  ## Row and column names that disagree make isSymmetric() fail on attributes
  ## alone (the GTRE report of Round 15). Those names are dropped...
  N <- matrix(c(2, 1, 1, 4), 2, dimnames = list(c("a", "b"), c("c", "d")))
  expect_null(dimnames(.safe_symmetrize(N)))
  ## ...but names that agree are kept.
  K <- matrix(c(2, 1, 1, 4), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_identical(dimnames(.safe_symmetrize(K)), dimnames(K))
})

test_that(".format_formula() pads a pipe formula to three parts with 1", {
  expect_identical(format(.format_formula(y ~ x1 + x2)), "y ~ x1 + x2 | 1 | 1")
  expect_identical(format(.format_formula(y ~ x1 | z)), "y ~ x1 | z | 1")
  expect_identical(format(.format_formula(y ~ x1 | z | w)), "y ~ x1 | z | w")
})

test_that(".parse_pipe_formula() splits the parts and names the response", {
  p1 <- .parse_pipe_formula(y ~ x1 + x2)
  expect_identical(p1$n_parts, 1L)
  expect_identical(p1$y_var, "y")
  expect_null(p1$formula_z)
  expect_null(p1$formula_zp)
  expect_identical(all.vars(p1$formula_x), c("y", "x1", "x2"))

  p3 <- .parse_pipe_formula(log(out) ~ x1 | z1 + z2 | w)
  expect_identical(p3$n_parts, 3L)
  expect_identical(p3$y_var, "out")
  expect_identical(all.vars(p3$formula_z), c("out", "z1", "z2"))
  expect_identical(all.vars(p3$formula_zp), c("out", "w"))
})
