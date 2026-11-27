## Arguments from sfaR and frontier are rejected with the sfa spelling
## (R/foreign_args.R). Every case errors before any fitting, so these are fast.

.fa_data <- function() {
  set.seed(1)
  data.frame(y = rnorm(20), x1 = rnorm(20), x2 = rnorm(20), id = rep(1:5, each = 4))
}

test_that("sfaR's udist and S translate to model_name and inefdec", {
  d <- .fa_data()
  expect_error(sfm(y ~ x1, data = d, udist = "hnormal"),
               'sfaR syntax. In sfa write model_name = "NHN"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, udist = "tnormal"), 'model_name = "NTN"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, udist = "exponential"), 'model_name = "NE"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, S = -1), "inefdec = FALSE", fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, S = 1), "inefdec = TRUE", fixed = TRUE)
  ## The value is looked up in the caller's frame.
  u <- "weibull"
  expect_error(sfm(y ~ x1, data = d, udist = u), 'model_name = "NW"', fixed = TRUE)
  ## An unevaluable value still gets the general advice.
  expect_error(sfm(y ~ x1, data = d, udist = no_such_object_xyz),
               "sfa sets the distribution through model_name", fixed = TRUE)
})

test_that("frontier's ineffDecrease, truncNorm and timeEffect translate", {
  d <- .fa_data()
  expect_error(sfm(y ~ x1, data = d, ineffDecrease = FALSE), "inefdec = FALSE", fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, truncNorm = TRUE), 'model_name = "NTN"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, timeEffect = TRUE), 'psfm(..., model_name = "BC92"', fixed = TRUE)
  expect_error(psfm(y ~ x1, data = d, individual = "id", timeEffect = TRUE),
               'model_name = "BC92"', fixed = TRUE)
  expect_error(psfm(y ~ x1, data = d, individual = "id", timeEffect = FALSE),
               'model_name = "PL80"', fixed = TRUE)
  expect_error(psfm(y ~ x1, data = d, individual = "id", truncNorm = TRUE),
               "no truncated-normal version", fixed = TRUE)
})

test_that("renames carry the value over, and name the right entry point", {
  d <- .fa_data()
  expect_error(lcsfm(y ~ x1, data = d, lcmClasses = 3), "In sfa write n_class = 3.", fixed = TRUE)
  expect_error(lcsfm(y ~ x1, data = d, thet = ~ x2), "In sfa write classhet = ~x2.", fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, lcmClasses = 3), "lcsfm() takes it as n_class", fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, simType = "sobol"), 'sim_type = "sobol"', fixed = TRUE)
  expect_error(zsfm(y ~ x1, data = d, seed = 1), "zsfm() has no equivalent", fixed = TRUE)
  expect_error(zsfm(y ~ x1, data = d, printInfo = TRUE), "verbose = TRUE", fixed = TRUE)
})

test_that("tuning arguments with no equivalent say so; unknown ones stay unused-argument errors", {
  d <- .fa_data()
  expect_error(sfm(y ~ x1, data = d, hessianType = 2), "no sfa equivalent; remove it", fixed = TRUE)
  expect_error(psfm(y ~ x1, data = d, individual = "id", gridSize = 0.1),
               "frontier syntax with no sfa equivalent", fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, foo = 3), "^unused argument \\(foo = 3\\)$")
  expect_error(sfm(y ~ x1, data = d, foo = 3, bar = 4), "unused arguments (foo = 3, bar = 4)", fixed = TRUE)
  ## Mixed: both the translation and the unknown argument are reported.
  e <- tryCatch(sfm(y ~ x1, data = d, udist = "hnormal", foo = 3), error = identity)
  expect_match(conditionMessage(e), 'model_name = "NHN"', fixed = TRUE)
  expect_match(conditionMessage(e), "unused argument (foo = 3)", fixed = TRUE)
  ## The error is attributed to the entry point's call, as R's own is.
  expect_identical(conditionCall(e)[[1]], as.name("sfm"))
})

test_that("partial matching of sfa's own arguments is unaffected by `...`", {
  d <- .fa_data()
  ## `model =` partially matches model_name; the call reaches psfm()'s own check.
  expect_error(psfm(y ~ x1, data = d, model = "TRE"), "needs `individual`", fixed = TRUE)
})

test_that("distribution names given as model_name point at the sfm model", {
  d <- .fa_data()
  expect_error(sfm(y ~ x1, data = d, model_name = "hnormal"),
               'names a distribution, not an sfa model. In sfa write model_name = "NHN"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, model_name = "Half-Normal"), 'model_name = "NHN"', fixed = TRUE)
  expect_error(sfm(y ~ x1, data = d, model_name = "truncated normal"), 'model_name = "NTN"', fixed = TRUE)
  ## Short forms keep the edit-distance suggestion they had before.
  expect_error(sfm(y ~ x1, data = d, model_name = "HN"), 'Did you mean "NHN"', fixed = TRUE)
})

test_that("the translation table only points at arguments and models that exist", {
  ## Every distribution alias names an sfm() model with a label.
  m <- unique(unname(sfa:::.SFA_DIST_ALIASES))
  expect_true(all(m %in% eval(formals(sfm)$model_name)))
  expect_true(all(m %in% names(sfa:::.SFA_DIST_LABEL)))
  ## The sfa arguments the advice names exist where the advice says they do.
  for (f in c("sfm", "psfm", "zsfm", "lcsfm")) {
    fm <- names(formals(get(f, envir = asNamespace("sfa"))))
    expect_true(all(c("inefdec", "start_val", "Method", "verbose", "maxit.optim", "...") %in% fm), info = f)
  }
  expect_true(all(c("sim_type", "sim_prime", "sim_burn", "sim_seed") %in% names(formals(sfm))))
  expect_true(all(c("n_class", "classhet") %in% names(formals(lcsfm))))
  expect_true(all(c("individual", "time") %in% names(formals(psfm))))
  expect_true(all(c("PL80", "BC92") %in% eval(formals(psfm)$model_name)))
  expect_true("logDepVar" %in% names(formals(efficiency)))
  ## No foreign name collides with an sfa argument (it would never reach `...`).
  foreign <- c(names(sfa:::.SFA_FOREIGN_ARGS), unlist(sfa:::.SFA_FOREIGN_NOEQUIV))
  for (f in c("sfm", "psfm", "zsfm", "lcsfm")) {
    fm <- names(formals(get(f, envir = asNamespace("sfa"))))
    expect_length(intersect(foreign, fm), 0)
  }
})
