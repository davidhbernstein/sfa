## NLN/NW: the composed density by quadrature rather than simulation (#65).

## Reference: adaptive integrate() in t = log u, where both inefficiency
## densities are smooth, split around the integrand's own mode and shifted by
## its maximum so a density far below double.xmin still compares.
nsml_ref <- function(e, sv, model, su, shp) {
  ld <- if (model == "NW") {
    function(u) dweibull(u, shape = shp, scale = su, log = TRUE)
  } else {
    function(u) dlnorm(u, meanlog = shp, sdlog = su, log = TRUE)
  }
  vapply(e, function(ei) {
    lg <- function(t) {
      v <- dnorm(ei + exp(t), 0, sv, log = TRUE) + ld(exp(t)) + t
      v[is.na(v)] <- -Inf
      v
    }
    tt <- seq(-120, 6, length.out = 20001)
    j <- which.max(lg(tt))
    om <- optimize(lg, tt[c(max(1, j - 1), min(length(tt), j + 1))],
      maximum = TRUE, tol = 1e-12)
    M <- om$objective
    br <- om$maximum + c(-60, -30, -15, -8, -4, -2, -1, -.5, -.25, -.1, -.03, 0,
      .03, .1, .25, .5, 1, 2, 4, 8)
    br <- sort(unique(c(br, seq(-120, 6, by = 3))))
    M + log(sum(vapply(seq_len(length(br) - 1L), function(i) {
      integrate(function(t) exp(lg(t) - M), br[i], br[i + 1L],
        rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L,
        stop.on.error = FALSE)$value
    }, 0)))
  }, 0)
}

nsml_cases <- list(
  ## model, sigma_v, sigma_u, shape
  list("NW", 0.3, 0.3, 1.5),   # issue #65's design
  list("NW", 0.3, 3, 0.3),     # k < 1: f_U singular at u = 0
  list("NW", 0.02, 0.1, 5),    # narrow noise, peaked u
  list("NW", 1, 0.01, 0.6),
  list("NLN", 0.3, 0.5, -1),
  list("NLN", 0.02, 0.01, 0.5), # narrow noise AND a 1% lognormal
  list("NLN", 0.05, 0.1, 0)
)

test_that("the density matches adaptive integration to near machine precision", {
  es <- c(-3, -1.5, -0.6, -0.2, 0, 0.1, 0.5, 1.5)
  for (cs in nsml_cases) {
    q <- sfa:::.log_d_nsml(es, cs[[2]], cs[[1]], cs[[3]], cs[[4]])
    r <- nsml_ref(es, cs[[2]], cs[[1]], cs[[3]], cs[[4]])
    keep <- r > -50
    expect_equal(q[keep], r[keep], tolerance = 1e-9,
      label = paste(unlist(cs), collapse = " "))
  }
})

test_that("the density integrates to one over the residual", {
  ## A check that uses no reference integral of the same kind.
  for (cs in nsml_cases[c(1, 2, 5)]) {
    f <- function(e) exp(sfa:::.log_d_nsml(e, cs[[2]], cs[[1]], cs[[3]], cs[[4]]))
    tot <- integrate(f, -Inf, Inf, rel.tol = 1e-10, subdivisions = 1000L)$value
    expect_equal(tot, 1, tolerance = 1e-8, label = paste(unlist(cs), collapse = " "))
  }
})

test_that("the efficiency predictor is the posterior mean over the same nodes", {
  sv <- 0.3; su <- 0.3; k <- 1.5
  es <- c(-1, -0.3, 0, 0.4)
  eff <- sfa:::.nsml_eff(es, sv, "NW", su, k)
  ref <- vapply(es, function(ei) {
    g <- function(u) dnorm(ei + u, 0, sv) * dweibull(u, k, su)
    integrate(function(u) exp(-u) * g(u), 0, Inf, rel.tol = 1e-12)$value /
      integrate(g, 0, Inf, rel.tol = 1e-12)$value
  }, 0)
  refu <- vapply(es, function(ei) {
    g <- function(u) dnorm(ei + u, 0, sv) * dweibull(u, k, su)
    integrate(function(u) u * g(u), 0, Inf, rel.tol = 1e-12)$value /
      integrate(g, 0, Inf, rel.tol = 1e-12)$value
  }, 0)
  expect_equal(eff$exp_u_hat, ref, tolerance = 1e-9)
  expect_equal(eff$u_hat, refu, tolerance = 1e-9)
  ## Jensen: E[exp(-u)] >= exp(-E[u]).
  expect_true(all(eff$exp_u_hat >= exp(-eff$u_hat) - 1e-12))
})

test_that("extreme parameters give finite values and no warnings", {
  es <- c(-200, -10, -1, 0, 1, 10, 50)
  for (model in c("NW", "NLN")) {
    for (sv in c(1e-3, 10)) for (su in c(1e-3, 5)) {
      for (shp in if (model == "NW") c(0.1, 20) else c(-5, 2)) {
        expect_silent(v <- sfa:::.log_d_nsml(es, sv, model, su, shp))
        expect_false(anyNA(v))
        expect_true(all(v < Inf))
      }
    }
  }
})

test_that("sfm() integrates by quadrature by default and records it", {
  d <- data_gen_cs(N = 60, rand = 7, sig_u = 0.3, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5, k_w = 1.5)
  expect_error(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d,
    integration = "monte"), "should be one of")
  skip_on_cran()
  expect_warning(
    f <- sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 800),
    "`Nsim` is ignored"
  )
  expect_identical(f$integration, "quadrature")
  ## logLik() is the quadrature density summed at the estimates: there is no
  ## simulation error left between them.
  p <- coef(f)
  e <- d$y_pcs_wb - cbind(1, d$x1, d$x2) %*% p[c("(Intercept)", "x1", "x2")]
  expect_equal(as.numeric(logLik(f)),
    sum(sfa:::.log_d_nsml(e, p[["sigv"]], "NW", p[["sigu"]], p[["k"]])),
    tolerance = 1e-10)
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
})
