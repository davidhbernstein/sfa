## A108: a cost fit of -y is the mirror of the production fit of y. With the
## starting values oriented and the simulated firm-effect draws following the
## orientation, the cost objective at mirrored parameters IS the production
## objective, and the cost starts are the mirrored production starts, so the
## whole procedure is equivariant: this asserts that property of the
## procedure, not where any optimizer lands.

.mirror <- function(model, formula, d, y, ...) {
  dn <- d; dn[[y]] <- -dn[[y]]
  p <- suppressWarnings(psfm(formula, model_name = model, data = d,
    individual = "name", ...))
  cst <- suppressWarnings(psfm(formula, model_name = model, data = dn,
    individual = "name", inefdec = FALSE, ...))
  list(p = p, c = cst)
}

test_that("TRE, TRE_Z and GTRE_Z cost fits mirror the production fits (A108)", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 5, N = 30, rand = 4, sig_u = 1, sig_v = 0.3,
    sig_r = 0.3, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  cases <- list(
    list("TRE", y_tre ~ x1 + x2, "y_tre", fr = c("(Intercept)", "x1", "x2")),
    list("TRE_Z", y_tre ~ x1 + x2 | z_gtre, "y_tre", fr = c("(Intercept x)", "x1", "x2")),
    list("GTRE_Z", y_gtre_z ~ x1 + x2 | z_gtre | zp_gtre, "y_gtre_z", fr = c("(Intercept x)", "x1", "x2")))
  for (k in cases) {
    ## A108 is about the simulated draws, so pin simulation (TRE/TRE_Z default
    ## to quadrature since 1.3.0, where the mirror is exact by construction).
    m <- .mirror(k[[1]], k[[2]], d, k[[3]], halton_num = 30, integration = "simulate")
    bp <- m$p$out[, "par"]; bc <- m$c$out[, "par"]
    expect_equal(unname(bc[k$fr]), unname(-bp[k$fr]), tolerance = 1e-5, info = k[[1]])
    sc <- setdiff(names(bp), k$fr)
    ## 1e-4 for the scales: TRE_Z's sigr sits on its zero boundary here, and the
    ## flat direction there amplifies last-bit differences between y and -y
    ## arithmetic to ~3e-5 in the z block. Before A108 these were 1e-2 to 1e-1.
    expect_equal(unname(abs(bc[sc])), unname(abs(bp[sc])), tolerance = 1e-4, info = k[[1]])
    expect_equal(as.numeric(logLik(m$c)), as.numeric(logLik(m$p)), tolerance = 1e-7, info = k[[1]])
  }
})
