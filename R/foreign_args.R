## Arguments that other R frontier packages (sfaR, frontier) use, and what sfa
## calls them. Entry points take `...` only so these can be caught and answered
## with the sfa spelling; anything else in `...` is still an unused-argument
## error.

## Distribution names, as sfaR's `udist` and common prose spell them, mapped to
## sfm() model names. Keys are normalised by .sfa_norm_dist().
.SFA_DIST_ALIASES <- c(
  hnormal = "NHN", halfnormal = "NHN",
  exponential = "NE", exp = "NE",
  tnormal = "NTN", truncnormal = "NTN", truncatednormal = "NTN",
  rayleigh = "NR", uniform = "NU", gamma = "NG", lognormal = "NLN",
  weibull = "NW", genexponential = "NGE", generalizedexponential = "NGE",
  tslaplace = "TSL", truncatedskewlaplace = "TSL"
)

.SFA_DIST_LABEL <- c(
  NHN = "half-normal", NE = "exponential", NTN = "truncated normal",
  NR = "Rayleigh", NU = "uniform", NG = "gamma", NLN = "lognormal",
  NW = "Weibull", NGE = "generalized exponential", TSL = "truncated skew-Laplace"
)

.sfa_norm_dist <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

## The sfm() model name for a distribution alias, or NA.
.sfa_dist_alias <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    return(NA_character_)
  }
  unname(.SFA_DIST_ALIASES[.sfa_norm_dist(x)])
}

.sfa_dist_advice <- function(v, fun) {
  m <- .sfa_dist_alias(v)
  if (fun == "sfm" && !is.na(m)) {
    sprintf('In sfa write model_name = "%s" (normal noise, %s inefficiency).',
            m, .SFA_DIST_LABEL[[m]])
  } else {
    sprintf("sfa sets the distribution through model_name; see ?%s.", fun)
  }
}

## One advice string per foreign argument. `v` is the argument's value, or
## NULL when it could not be evaluated; `fun` is the sfa entry point; `has`
## says whether that entry point has a given argument.
.SFA_FOREIGN_ARGS <- list(
  udist = list(pkg = "sfaR", f = function(v, fun, has) .sfa_dist_advice(v, fun)),
  dist = list(pkg = NA, f = function(v, fun, has) .sfa_dist_advice(v, fun)),
  distribution = list(pkg = NA, f = function(v, fun, has) .sfa_dist_advice(v, fun)),
  S = list(pkg = "sfaR", f = function(v, fun, has) {
    if (identical(as.numeric(v), -1)) {
      "In sfa write inefdec = FALSE for a cost frontier."
    } else if (identical(as.numeric(v), 1)) {
      "In sfa write inefdec = TRUE for a production frontier (the default)."
    } else {
      "In sfa write inefdec = TRUE for a production frontier (sfaR's S = 1) or inefdec = FALSE for a cost frontier (S = -1)."
    }
  }),
  ineffDecrease = list(pkg = "frontier", f = function(v, fun, has) {
    if (isFALSE(v)) {
      "In sfa write inefdec = FALSE for a cost frontier."
    } else {
      "In sfa write inefdec = TRUE (the default) for a production frontier, FALSE for a cost one."
    }
  }),
  truncNorm = list(pkg = "frontier", f = function(v, fun, has) {
    if (fun == "psfm") {
      if (isTRUE(v)) {
        'sfa\'s "PL80" and "BC92" are half-normal; there is no truncated-normal version of them.'
      } else {
        "Remove it: sfa's panel error-components models are half-normal."
      }
    } else if (fun == "sfm") {
      if (isTRUE(v)) 'In sfa write model_name = "NTN".' else 'Remove it: model_name = "NHN" is the default.'
    } else {
      sprintf("sfa sets the distribution through model_name; see ?%s.", fun)
    }
  }),
  timeEffect = list(pkg = "frontier", f = function(v, fun, has) {
    m <- if (isTRUE(v)) "BC92" else "PL80"
    if (fun == "psfm") {
      sprintf('In sfa write model_name = "%s", with individual = and time = naming the firm and period columns.', m)
    } else {
      sprintf('That is a panel model: use psfm(..., model_name = "%s", individual = "<firm column>", time = "<period column>").', m)
    }
  }),
  startVal = list(pkg = "frontier", f = function(v, fun, has) "In sfa write start_val ="),
  logDepVar = list(pkg = "sfaR", f = function(v, fun, has) {
    "The fit does not use it: pass it to efficiency(fit, logDepVar = ) instead."
  }),
  method = list(pkg = "sfaR", f = function(v, fun, has) {
    "sfa's final optim() stage takes Method = (default \"L-BFGS-B\"); sfaR's method names do not carry over."
  }),
  itermax = list(pkg = "sfaR", f = function(v, fun, has) {
    "In sfa set maxit.optim =, and maxit.bobyqa = / maxit.psoptim = for the earlier stages."
  }),
  printInfo = list(pkg = "sfaR", f = function(v, fun, has) "In sfa write verbose = TRUE."),
  printIter = list(pkg = "frontier", f = function(v, fun, has) "In sfa write verbose = TRUE."),
  simType = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("sim_type", fun, has)),
  prime = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("sim_prime", fun, has)),
  burn = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("sim_burn", fun, has)),
  seed = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("sim_seed", fun, has)),
  lcmClasses = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("n_class", fun, has, "lcsfm")),
  thet = list(pkg = "sfaR", f = function(v, fun, has) .sfa_rename("classhet", fun, has, "lcsfm")),
  subset = list(pkg = NA, f = function(v, fun, has) {
    "Subset the data before the call instead: data = subset(d, <condition>)."
  })
)

## Tuning arguments with no sfa counterpart.
.SFA_FOREIGN_NOEQUIV <- list(
  sfaR = c("hessianType", "tol", "gradtol", "stepmax", "qac", "whichStart",
           "initAlg", "initIter", "randStart"),
  frontier = c("muBound", "bignum", "searchStep", "searchTol", "searchScale",
               "gridSize", "gridDouble", "restartMax", "restartFactor")
)

.sfa_rename <- function(to, fun, has, where = NULL) {
  if (has(to)) {
    sprintf("In sfa write %s =", to)
  } else if (!is.null(where)) {
    sprintf("%s() has no such argument; %s() takes it as %s =.", fun, where, to)
  } else {
    sprintf("%s() has no equivalent.", fun)
  }
}

## Called first thing in an entry point as
##   .sfa_reject_dots(match.call(expand.dots = FALSE)$..., "sfm", parent.frame())
## Errors if `...` holds anything, translating what it recognises.
.sfa_reject_dots <- function(dots, fun, env) {
  if (!length(dots)) {
    return(invisible(NULL))
  }
  call <- sys.call(-1L)
  fmls <- names(formals(get(fun, envir = asNamespace("sfa"))))
  has <- function(a) a %in% fmls
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))

  show <- function(i) {
    e <- paste(deparse(dots[[i]], width.cutoff = 60L)[1L], collapse = "")
    if (nzchar(nms[i])) paste(nms[i], "=", e) else e
  }
  advice <- character(0)
  unknown <- character(0)
  for (i in seq_along(dots)) {
    nm <- nms[i]
    hit <- if (nzchar(nm)) .SFA_FOREIGN_ARGS[[nm]] else NULL
    noeq <- if (nzchar(nm)) names(Filter(function(a) nm %in% a, .SFA_FOREIGN_NOEQUIV)) else character(0)
    if (!is.null(hit)) {
      v <- tryCatch(eval(dots[[i]], env), error = function(e) NULL)
      who <- if (is.na(hit$pkg)) "is not an sfa argument." else paste0("is ", hit$pkg, " syntax.")
      a <- hit$f(v, fun, has)
      ## A bare rename ("In sfa write n_class =") carries the user's value over.
      if (endsWith(a, " =")) {
        a <- paste0(a, " ", paste(deparse(dots[[i]], width.cutoff = 60L)[1L], collapse = ""), ".")
      }
      advice <- c(advice, paste0("`", show(i), "` ", who, " ", a))
    } else if (length(noeq)) {
      advice <- c(advice, paste0("`", show(i), "` is ", noeq[1L],
                                 " syntax with no sfa equivalent; remove it."))
    } else {
      unknown <- c(unknown, show(i))
    }
  }
  msg <- character(0)
  if (length(advice)) {
    msg <- c(paste0(fun, "() does not take ",
                    if (length(advice) == 1L) "this argument" else "these arguments",
                    " from another package:"),
             paste0("  ", advice))
  }
  if (length(unknown)) {
    msg <- c(msg, paste0(if (length(advice)) "  ", "unused argument", if (length(unknown) > 1L) "s", " (",
                         paste(unknown, collapse = ", "), ")"))
  }
  stop(simpleError(paste(msg, collapse = "\n"), call = call))
}
