## tt_boundary_report(): is a collapsed sigma_v the maximum, or did the
## optimizer stop short?  Gap A20.
##
## Chris Parmeter asked for the exercise that settled ZISF's boundary question
## to be repeated for the two-tier models: refit from several starts, score
## every solution on ONE likelihood, and see whether the collapsed point is
## actually the best. That is all this does.
##
## IT IS NOT A LIKELIHOOD-RATIO TEST and nothing here reports a p-value.
## Chris Parmeter, 2026-09-11: "I am fine with a likelihood comparison as a
## diagnostic, but I would avoid describing it as a conventional LR test
## unless we establish the appropriate null distribution given the boundary
## problem." The null sits ON the boundary of the parameter space, so the
## usual chi-square reference does not apply, and no appropriate null
## distribution has been established for this model. The output is a table of
## log-likelihoods and a statement about which is largest.

tt_boundary_report <- function(object, data,
                               sigma_v_starts = c(0.05, 0.3, 1),
                               tol = 1e-3, tol_loglik = 1e-4, ...) {
  if (!inherits(object, "sfareg") ||
      !identical(object$model_name, "TTHN") &&
      !identical(object$model_name, "TTNE")) {
    stop("`object` must be a ttsfm() fit with model_name \"TTHN\" or ",
      "\"TTNE\". The diagnostic is about the two-tier noise scale and has no ",
      "meaning elsewhere.",
      call. = FALSE
    )
  }
  if (missing(data) || is.null(data)) {
    stop("`data` must be supplied: ttsfm() does not retain the data it was ",
      "fitted to, and every refit here needs it.",
      call. = FALSE
    )
  }
  if (!is.numeric(sigma_v_starts) || !length(sigma_v_starts) ||
      any(!is.finite(sigma_v_starts)) || any(sigma_v_starts <= 0)) {
    stop("`sigma_v_starts` must be a vector of positive, finite starting ",
      "values for sigma_v.",
      call. = FALSE
    )
  }

  p0 <- object$out[, "par"]
  nm <- names(p0)
  i_v <- match("sigma_v", nm)
  if (is.na(i_v)) {
    stop("This fit does not report a `sigma_v` row. Fits made before 1.2.1 ",
      "labelled it `sigv` and carried it on the LOG scale; refit under the ",
      "current version.",
      call. = FALSE
    )
  }

  ## Back to the raw optimizer scale ttsfm()'s `start_val` expects: log
  ## sigma_v, then the two determinant blocks as they already are.
  raw <- p0
  raw[i_v] <- log(pmax(p0[i_v], .Machine$double.eps))
  i_u <- match("sigma_u", nm)
  i_w <- match("sigma_w", nm)
  z_link <- if (is.null(object$z_link)) "sd" else object$z_link
  to_eta <- function(s) if (identical(z_link, "sd")) log(s) else 2 * log(s)
  if (!is.na(i_u)) raw[i_u] <- to_eta(pmax(p0[i_u], .Machine$double.eps))
  if (!is.na(i_w)) raw[i_w] <- to_eta(pmax(p0[i_w], .Machine$double.eps))

  ll0 <- suppressWarnings(as.numeric(stats::logLik(object)))
  rows <- data.frame(start = "as returned", sigma_v_start = NA_real_,
    sigma_v = unname(p0[i_v]), logLik = ll0, stringsAsFactors = FALSE)

  for (sv in sigma_v_starts) {
    st <- raw
    st[i_v] <- log(sv)
    fit <- tryCatch(
      suppressWarnings(ttsfm(object$formula, data = data,
        model_name = object$model_name, start_val = unname(st), ...
      )),
      error = function(e) NULL
    )
    rows <- rbind(rows, data.frame(
      start = sprintf("sigma_v = %g", sv), sigma_v_start = sv,
      sigma_v = if (is.null(fit)) NA_real_ else unname(fit$out[["sigma_v", "par"]]),
      logLik = if (is.null(fit)) NA_real_ else {
        suppressWarnings(as.numeric(stats::logLik(fit)))
      },
      stringsAsFactors = FALSE
    ))
  }

  ## The ordinary normal linear model on the same data, as the floor: a
  ## collapsed sigma_v is NOT the same thing as "there is no frontier here",
  ## and the gap to this row is what says so.
  fx <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
  lmfit <- tryCatch(stats::lm(fx, data = data), error = function(e) NULL)
  ll_lm <- if (is.null(lmfit)) NA_real_ else as.numeric(stats::logLik(lmfit))

  ## Ties go to the returned fit. Restarts that land on the same optimum
  ## differ in the last few digits, and which.max() would then report a
  ## "higher" likelihood of 0 log units -- which is how this printed on its
  ## first outing. `tol_loglik` is what counts as meaningfully higher.
  lv <- replace(rows$logLik, is.na(rows$logLik), -Inf)
  best <- if (max(lv) - lv[1] > tol_loglik) which.max(lv) else 1L
  collapsed <- is.finite(rows$sigma_v[1]) && rows$sigma_v[1] < tol

  res <- list(
    table = rows, best = rows$start[best],
    returned_is_best = isTRUE(best == 1L),
    gap = if (length(rows$logLik) > 1) {
      max(rows$logLik[-1], na.rm = TRUE) - ll0
    } else NA_real_,
    tol_loglik = tol_loglik,
    collapsed = collapsed, tol = tol,
    lm_logLik = ll_lm, lm_gap = ll0 - ll_lm,
    model_name = object$model_name, call = match.call()
  )
  class(res) <- "sfa_tt_boundary"
  res
}

print.sfa_tt_boundary <- function(x, ...) {
  cat("Two-tier noise-scale boundary report (", x$model_name, ")\n", sep = "")
  cat("A LIKELIHOOD COMPARISON, NOT A TEST -- no p-value is reported, and none\n")
  cat("should be inferred: the null lies on the boundary of the parameter space.\n\n")
  tb <- x$table
  tb$logLik <- round(tb$logLik, 4)
  tb$sigma_v <- signif(tb$sigma_v, 4)
  print(tb[, c("start", "sigma_v", "logLik")], row.names = FALSE)
  cat("\n")
  if (x$collapsed) {
    cat("The returned fit has sigma_v <", format(x$tol),
      "-- it has collapsed.\n", sep = " ")
  }
  if (isTRUE(x$returned_is_best)) {
    cat("No restart found a higher likelihood.")
    if (x$collapsed) {
      cat(" On this evidence the collapse IS the maximum,\nnot a failure of the optimizer.")
    }
    cat("\n")
  } else {
    cat("A restart found a HIGHER likelihood, by",
      format(round(x$gap, 4)), "log units:", x$best, "\n")
    cat("The returned fit is not the maximum. Refit from that start.\n")
  }
  if (is.finite(x$lm_gap)) {
    cat("\nAgainst the ordinary normal linear model on the same data: the fit\n",
      "is ", format(round(x$lm_gap, 2)), " log units ",
      if (x$lm_gap > 0) "better" else "worse",
      " (lm logLik ", format(round(x$lm_logLik, 2)), ").",
      if (x$lm_gap > 0 && isTRUE(x$collapsed)) {
        " A collapsed sigma_v is not\nthe same thing as no frontier."
      } else "",
      "\n", sep = ""
    )
  }
  invisible(x)
}
