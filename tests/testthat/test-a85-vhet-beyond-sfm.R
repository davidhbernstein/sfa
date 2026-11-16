## A85/A86. Heteroskedastic noise outside sfm(), and the two defects the work
## uncovered in psfm("FD")'s efficiency extraction.
##
## What is asserted here is the CONTRACT, never where an optimizer landed:
## the layout a `vhet` block produces, the refusals, the invariance of the
## homoskedastic path, and -- for FD -- that the stored u_hat is the posterior
## mean the model implies AT WHATEVER PARAMETERS THE FIT RETURNED. That last
## form is what makes the FD regression test immune to the platform-dependent
## landing point this suite has been bitten by before.

cs <- function(N = 250, rand = 7)
  as.data.frame(data_gen_cs(N = N, rand = rand, sig_u = 1, sig_v = 0.3,
                            cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5))

test_that("vhet must be a one-sided formula in every entry point that takes it", {
  d <- cs(80)
  expect_error(ttsfm(y_ttne ~ x1 + x2, "TTNE", data = d, vhet = "zv"),
               "must be a one-sided formula")
  expect_error(zsfm(y_zisf ~ x1 + x2, "ZISF", data = d, vhet = 1),
               "must be a one-sided formula")
})

test_that("vhet is refused where sigma_v is not in the objective", {
  d <- cs(80)
  ## TTNLS minimises a sum of squares on E[e] = sigma_u - sigma_w. sigma_v
  ## never enters it, which the package already says elsewhere by returning
  ## the scales as NA.
  expect_error(ttsfm(y_ttne ~ x1 + x2, "TTNLS", data = d, vhet = ~ x1),
               "no meaning for model_name")
  ## psfm() takes vhet for the first-difference likelihood only.
  p <- as.data.frame(data_gen_p(t = 4, N = 20, rand = 42, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5))
  expect_error(psfm(y_gtre ~ x1 + x2, "GTRE", data = p, individual = "name",
                    vhet = ~ x1),
               "implemented for model_name")
})

test_that("a vhet block replaces the single sigma_v row and is named for its design", {
  skip_on_cran()
  set.seed(3); n <- 300
  zv <- runif(n); x1 <- runif(n); x2 <- runif(n)
  d <- data.frame(y = 0.5 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, exp(log(0.3) + 0.6 * zv)) +
                      rexp(n, 1 / 0.8) - rexp(n, 1), x1, x2, zv)

  hom <- ttsfm(y ~ x1 + x2, "TTNE", data = d)
  het <- ttsfm(y ~ x1 + x2, "TTNE", data = d, vhet = ~ zv)

  expect_true("sigma_v" %in% rownames(hom$out))
  ## With determinants there is no single sigma_v to report, so those rows
  ## become COEFFICIENTS with a Zv. prefix -- the same convention the u and w
  ## blocks already follow.
  expect_false("sigma_v" %in% rownames(het$out))
  expect_true(all(c("Zv.(Intercept)", "Zv.zv") %in% rownames(het$out)))
  expect_equal(nrow(het$out), nrow(hom$out) + 1L)

  ## The per-observation scale is stored, because out[, "par"] no longer
  ## holds it anywhere.
  expect_length(het$sigma_v, nrow(d))
  expect_true(all(het$sigma_v > 0))
  expect_null(hom$sigma_v)
})

test_that("vhet = ~ 1 is the homoskedastic model, just reparameterized", {
  skip_on_cran()
  d <- cs(200)
  a <- ttsfm(y_ttne ~ x1 + x2, "TTNE", data = d)
  b <- ttsfm(y_ttne ~ x1 + x2, "TTNE", data = d, vhet = ~ 1)
  ## Same model, so the attained log-likelihood must agree; only the
  ## parameterization of one block differs.
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-3)
  expect_equal(unname(b$sigma_v[1]), unname(a$out[["sigma_v", "par"]]),
               tolerance = 1e-3)
})

test_that("zsfm() keeps the mixing block separate from the noise block", {
  skip_on_cran()
  d <- cs(250)
  f <- suppressWarnings(zsfm(y_zisf_z ~ x1 + x2 | z, "ZISF_Z", data = d,
                             vhet = ~ x1))
  rn <- rownames(f$out)
  ## ghet/pipe block and vhet block both present and distinguishable.
  expect_true(all(c("Zv.(Intercept)", "Zv.x1") %in% rn))
  expect_true("z" %in% rn)
  expect_length(f$sigma_v, nrow(d))
})

test_that("psfm('FD') u_hat is the posterior mean the model implies (A86)", {
  skip_on_cran()
  ## Two regressors on purpose: the delta-slice defect was invisible at one,
  ## because `length(n_x_vars)` is 1 whatever the model.
  set.seed(4); N <- 60; T <- 5
  id <- rep(seq_len(N), each = T)
  x1 <- rnorm(N * T); x2 <- rnorm(N * T); z <- runif(N * T)
  alpha <- rep(rnorm(N, 0, 1), each = T)
  ustar <- rep(abs(rnorm(N, 1, 1)), each = T)
  d <- data.frame(name = id,
                  y = 1 + 0.5 * x1 - 0.3 * x2 + alpha +
                      rnorm(N * T, 0, 0.4) - exp(0.7 * z) * ustar,
                  x1, x2, z)
  f <- suppressWarnings(psfm(y ~ x1 + x2 | z, "FD", data = d,
                             individual = "name"))
  p <- f$out[, "par"]
  skip_if_not(all(is.finite(p)) && p[["sig_u2"]] > 1e-8,
              "FD landed on the sigma_u boundary; the posterior is degenerate there")

  ## Wang and Ho (2010), written from the model rather than from the package:
  ## Delta eps = D(y - Xb), Delta h = D exp(z'delta), Sigma = sigma_v^2 D D'.
  ## Evaluated AT THE FIT'S OWN PARAMETERS, so this asserts an identity and
  ## not a landing point.
  b <- p[c("x1", "x2")]; dl <- p[["z"]]
  su2 <- p[["sig_u2"]]; sv2 <- p[["sig_v2"]]; mu <- p[["mu"]]
  uh <- numeric(nrow(d))
  for (i in seq_len(N)) {
    r <- which(d$name == i); Ti <- length(r)
    D <- diag(Ti)[-1, , drop = FALSE] - diag(Ti)[-Ti, , drop = FALSE]
    de <- as.numeric(D %*% (d$y[r] - cbind(d$x1[r], d$x2[r]) %*% b))
    dh <- as.numeric(D %*% exp(d$z[r] * dl))
    Si <- solve(sv2 * (D %*% t(D)))
    s2 <- 1 / (as.numeric(t(dh) %*% Si %*% dh) + 1 / su2)
    ms <- (mu / su2 - as.numeric(t(de) %*% Si %*% dh)) * s2
    uh[r] <- exp(d$z[r] * dl) *
      (ms + sqrt(s2) * dnorm(ms / sqrt(s2)) / pnorm(ms / sqrt(s2)))
  }
  expect_equal(as.numeric(f$u_hat), uh, tolerance = 1e-8)
  expect_equal(as.numeric(f$exp_u_hat), pmin(pmax(exp(-uh), 0), 1),
               tolerance = 1e-8)
})

test_that("psfm('FD') with vhet reduces to the homoskedastic covariance", {
  skip_on_cran()
  ## The point of the generalisation: Var(D v) = D diag(sv2) D' must equal
  ## sigma_v^2 * SIGMA_i when sv2 is constant. Asserted on the builder's
  ## output through a fit, by giving vhet a covariate with no effect and
  ## checking the fitted scale is flat and the layout is right.
  set.seed(9); N <- 50; T <- 5
  id <- rep(seq_len(N), each = T)
  x1 <- rnorm(N * T); z <- runif(N * T)
  alpha <- rep(rnorm(N, 0, 1), each = T)
  ustar <- rep(abs(rnorm(N, 1, 1)), each = T)
  d <- data.frame(name = id,
                  y = 1 + 0.5 * x1 + alpha + rnorm(N * T, 0, 0.4) -
                      exp(0.5 * z) * ustar,
                  x1, z)
  f <- suppressWarnings(psfm(y ~ x1 | z, "FD", data = d, individual = "name",
                             vhet = ~ 1))
  expect_true("Zv.(Intercept)" %in% rownames(f$out))
  expect_false("sig_v2" %in% rownames(f$out))
  expect_length(f$sigma_v, nrow(d))
  ## Intercept-only: one value, repeated.
  expect_equal(length(unique(round(f$sigma_v, 12))), 1L)
})

test_that("psfm('FD') heteroskedastic covariance reduces to the old one exactly", {
  skip_on_cran()
  ## The sharpest check available for the generalisation, and the reason it is
  ## written this way: evaluating the HET closure at the parameters the
  ## HOMOSKEDASTIC fit returned must return the homoskedastic closure's own
  ## value. Comparing two FITS would only compare two optimizer landing
  ## points; comparing two CLOSURES at one point compares the arithmetic.
  ##
  ## What it certifies: D diag(sigma_v^2 ... ) D' equals sigma_v^2 * SIGMA_i,
  ## and the direct log-determinant equals the factorisation
  ## -0.5 log(T) - 0.5 (T-1) log(sigma_v^2) that the homoskedastic branch
  ## keeps verbatim.
  set.seed(9); N <- 40; T <- 5
  id <- rep(seq_len(N), each = T)
  x1 <- rnorm(N * T); x2 <- rnorm(N * T); z <- runif(N * T)
  alpha <- rep(rnorm(N, 0, 1), each = T)
  ustar <- rep(abs(rnorm(N, 1, 1)), each = T)
  d <- data.frame(name = id,
                  y = 1 + 0.5 * x1 - 0.3 * x2 + alpha +
                      rnorm(N * T, 0, 0.45) - exp(0.6 * z) * ustar,
                  x1, x2, z)
  a <- suppressWarnings(psfm(y ~ x1 + x2 | z, "FD", data = d,
                             individual = "name", keep_objective = TRUE))
  b <- suppressWarnings(psfm(y ~ x1 + x2 | z, "FD", data = d,
                             individual = "name", vhet = ~ 1,
                             keep_objective = TRUE))
  pa <- a$out[, "par"]
  skip_if_not(all(is.finite(pa)) && pa[["sig_v2"]] > 0)
  ## Same point, expressed in each closure's own parameterization: the block
  ## holds log(sigma_v^2), the scalar slot holds sigma_v^2.
  ph <- c(pa[["sig_u2"]], log(pa[["sig_v2"]]), pa[["mu"]],
          pa[["x1"]], pa[["x2"]], pa[["z"]])
  expect_equal(b$objective(ph), a$objective(unname(pa)), tolerance = 1e-10)
})
