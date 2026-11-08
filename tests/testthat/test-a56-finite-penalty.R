## A56. A draw the likelihood cannot evaluate was penalised by
## -sqrt(double.xmax / n) per observation, ~-7.7e152. NG and NNAK fits whose
## scale or shape sits on its lower bound have optim()'s Hessian stencil step
## across that bound, so the penalty entered the Hessian and every standard
## error was lost. The bound coordinate's SE was a function of the penalty
## (1.6e-6 at 1e4, 1.6e-7 at 1e6), so no constant fixes it: it gets none, and
## the rest come from the Hessian over the free coordinates.

test_that("the stencil rule flags only coordinates optim() would step across", {
  lo <- c(1e-7, 1e-7, -Inf, 0)
  expect_identical(sfa:::.hessian_stencil_at_bound(c(3.9e-6, 0.5, -10, 0.0021), lo),
    c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(sfa:::.hessian_stencil_at_bound(c(0.001, 1, 1, 0.0019), lo),
    c(TRUE, FALSE, FALSE, TRUE))
  ## A lower-bound vector of another length cannot be read: flag nothing.
  expect_false(any(sfa:::.hessian_stencil_at_bound(1:3, c(0, 0))))
})

test_that("no guard penalises with a power of double.xmax", {
  for (f in c("sfm", "lcsfm", "zsfm")) {
    src <- paste(deparse(get(f, envir = asNamespace("sfa"))), collapse = "\n")
    expect_false(grepl("double.xmax/length(like)", gsub(" ", "", src), fixed = TRUE),
      label = f)
  }
})

test_that("a coordinate on its bound loses only its own standard error", {
  skip_on_cran()
  for (cs in list(c("NG", "y_pcs_g"), c("NNAK", "y_pcs_nak"))) {
    for (r in 1:4) {
      d <- data_gen_cs(N = 150, rand = r, cons = .5, beta1 = .5, beta2 = .5,
        sig_u = 1, sig_v = .3, mu = .5, a = 5)
      w <- character(0)
      f <- withCallingHandlers(
        sfm(stats::as.formula(paste(cs[2], "~ x1 + x2")), data = d,
          model_name = cs[1]),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        })
      info <- paste(cs[1], "rand", r)
      se <- f$out[, "st_err"]
      if (length(f$se_at_bound)) {
        expect_true(all(is.na(se[f$se_at_bound])), info = info)
        expect_true(any(grepl("lies on its lower bound", w)), info = info)
        ## Every other SE is either finite or explained by $se_unresolved.
        rest <- setdiff(names(se), c(f$se_at_bound, f$se_unresolved$parameters))
        expect_true(all(is.finite(se[rest])), info = info)
      } else {
        expect_false(any(grepl("lies on its lower bound", w)), info = info)
      }
    }
  }
})
