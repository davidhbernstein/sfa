## Gap A36. .gtre_te() turned one firm's failed posterior inversion into stop(),
## so the whole GTRE_Z fit was lost -- it happened on the GitHub Actions Linux
## and Windows runners for a small boundary panel. That firm now gets NA
## efficiencies and a warning, and the estimates are returned.

test_that("GTRE_Z returns its fit with NA efficiencies for a firm whose posterior cannot be inverted", {
  skip_on_cran()
  ns <- asNamespace("sfa")
  orig <- get(".safe_linear_combo", envir = ns)
  on.exit(assignInNamespace(".safe_linear_combo", orig, ns = "sfa"), add = TRUE)
  calls <- 0L
  failing <- function(...) {
    calls <<- calls + 1L
    if (calls == 9L) stop("Unable to invert posterior system even after ridging.", call. = FALSE)
    orig(...)
  }
  assignInNamespace(".safe_linear_combo", failing, ns = "sfa")

  d <- as.data.frame(data_gen_p(t = 5, N = 60, rand = 11, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  w <- character(0)
  f <- withCallingHandlers(
    psfm(y_gtre_z ~ x1 + x2 | z_gtre | zp_gtre, model_name = "GTRE_Z", data = d,
      individual = "name", halton_num = 30, rand.gtre = 7),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  rows9 <- which(as.character(d$name) == as.character(unique(d$name)[9]))
  expect_true(all(is.finite(f$coefficients)))
  expect_true(all(is.na(f$U[rows9])))
  expect_true(all(is.finite(f$U[-rows9])))
  expect_true(any(grepl("firm\\(s\\) 9", w)))
})
