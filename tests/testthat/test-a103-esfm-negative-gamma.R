## A103: esfm(dist = "exponential") with gamma < 0. The documented density of u
## is the MIRROR, prop. to exp(-(B - u)/g) on [0, B], B = a0 |gamma|. Pinned
## against integrate() on the log scale, not against a fit.

.esfm_exp_neg <- function(g, sv) {
  a0 <- .HMS$exponential$a0; B <- a0 * g
  lf <- function(w) function(u) -(B - u) / g - log(g * (1 - exp(-a0))) +
    dnorm((w + u) / sv, log = TRUE) - log(sv)
  list(logdens = function(w) { f <- lf(w); m <- max(f(c(0, B, B / 2)))
         m + log(integrate(function(u) exp(f(u) - m), 0, B, rel.tol = 1e-13)$value) },
       te = function(w) { f <- lf(w); m <- max(f(seq(0, B, length.out = 200)))
         integrate(function(u) exp(-u + f(u) - m), 0, B, rel.tol = 1e-12)$value /
           integrate(function(u) exp(f(u) - m), 0, B, rel.tol = 1e-12)$value })
}

test_that("the gamma < 0 exponential log density is exact far into the right tail (A103b)", {
  ref <- .esfm_exp_neg(1.2, 0.9)
  w <- c(-8, -3, 0, 3, 7, 10, 20, 40)
  expect_equal(.esfm_logdens(w, -1.2, 0.9, "exponential"), vapply(w, ref$logdens, 1),
    tolerance = 1e-9)
  ## the old form gave -54.0 here, 200 log-units above the truth
  expect_lt(.esfm_logdens(20, -1.2, 0.9, "exponential"), -250)
})

test_that("the gamma < 0 exponential score integrates the mirrored posterior (A103a)", {
  ref <- .esfm_exp_neg(1.2, 0.9)
  w <- c(-8, -3, 0, 3, 7, 10, 20)
  fit <- list(gamma = -1.2, sigma_v = 0.9, dist = "exponential", inefdec = TRUE,
              out = cbind(par = c("(Intercept)" = 0)))
  X <- matrix(1, length(w), 1, dimnames = list(NULL, "(Intercept)"))
  expect_equal(.esfm_te(fit, w, X), vapply(w, ref$te, 1), tolerance = 1e-9)
  ## gamma > 0 is Battese-Coelli under N(-w - sv^2/g, sv^2) on [0, Inf)
  fit$gamma <- 1.2
  expect_equal(.esfm_te(fit, w, X), .te_battese_coelli(-w - 0.81 / 1.2, 0.9))
})

test_that("predict(type = 'efficiency') reads the same slot as efficiency()", {
  f <- structure(list(Efficiency = c(0.8, 0.6), model_name = "EXT-EXPONENTIAL"),
                 class = c("esfm", "sfareg"))
  expect_equal(predict(f, type = "efficiency"), c(0.8, 0.6))
  expect_equal(predict(f, type = "efficiency"), as.numeric(efficiency(f, type = "bc")))
})
