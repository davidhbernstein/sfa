## Gap A33. ttsfm("TTNE")'s E[exp(-u) | eps] (and so the M6 metric) used
## 0.5 * sig.v where the validated closed form (base code/ttsfm/2TierR.Rnw) has
## 0.5 * sig.v^2. The two agree only at sigma_v = 1. Checked here against the
## conditional expectation computed by numerical integration over (u, w) | eps.

.ttne_post_E <- function(lg, eps, sv, su, sw) {
  U <- 60 * max(su, sv) + abs(eps) + 5
  W <- 60 * max(sw, sv) + abs(eps) + 5
  lf <- function(u, w) dnorm((eps + u - w) / sv, log = TRUE) - log(sv) - u / su - log(su) - w / sw - log(sw)
  one <- function(g) {
    inner <- function(w) vapply(w, function(ww) integrate(function(u) exp(g(u, ww) + lf(u, ww)),
      0, U, rel.tol = 1e-10, subdivisions = 5000L)$value, numeric(1))
    integrate(inner, 0, W, rel.tol = 1e-9, subdivisions = 5000L)$value
  }
  one(lg) / one(function(u, w) 0)
}

test_that("TTNE conditional expectations match numerical integration", {
  skip_on_cran()
  set.seed(21)
  n <- 400; x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n, 0, 0.3) - rexp(n, 1 / 0.6) + rexp(n, 1 / 0.2)
  f <- suppressWarnings(ttsfm(y ~ x, model_name = "TTNE", data = data.frame(y, x)))
  ## Layout as metrics.ne() reads opt$par for a formula without pipes: frontier
  ## coefficients, then log sigma_v, then the u and w scale predictors, which
  ## the default z_link = "sd" maps through exp(); the residual is the raw
  ## y - x'b the likelihood uses.
  p <- f$opt$par
  sv <- exp(p[3]); su <- exp(p[4]); sw <- exp(p[5])
  expect_gt(abs(sv - 1), 0.3)
  ep <- y - p[1] - p[2] * x
  for (i in c(3, 50, 200)) {
    expect_equal(f$metrics$Eemu.cond[i], .ttne_post_E(function(u, w) -u, ep[i], sv, su, sw), tolerance = 1e-6)
    expect_equal(f$metrics$Eemw.cond[i], .ttne_post_E(function(u, w) -w, ep[i], sv, su, sw), tolerance = 1e-6)
  }
})
