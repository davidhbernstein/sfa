## Gap A32. FD's likelihood floored pnorm() at machine epsilon inside log() for
## two terms, so wherever mu*/sigma* or mu/sigma_u fell below about -8.1 the
## reported log-likelihood was not the likelihood: on this panel the old code
## reported -165.41 at its optimum where the exact value there is -166.10.

.fd_exact_loglik <- function(p, d) {
  su2 <- p[["sig_u2"]]; sv2 <- p[["sig_v2"]]; mu <- p[["mu"]]
  b <- p[["x_fd"]]; dl <- p[["z_fd"]]
  tot <- 0
  for (g in split(d, d$name)) {
    g <- g[order(g$year), ]
    Ti <- nrow(g)
    if (Ti < 2) next
    et <- diff(g$y_fd) - b * diff(g$x_fd)
    h <- diff(exp(dl * g$z_fd))
    S <- matrix(0, Ti - 1, Ti - 1); diag(S) <- 2
    if (Ti > 2) {
      S[cbind(1:(Ti - 2), 2:(Ti - 1))] <- -1
      S[cbind(2:(Ti - 1), 1:(Ti - 2))] <- -1
    }
    Si <- solve(sv2 * S)
    ss2 <- 1 / (as.numeric(t(h) %*% Si %*% h) + 1 / su2)
    ms <- (mu / su2 - as.numeric(t(et) %*% Si %*% h)) * ss2
    tot <- tot - 0.5 * (Ti - 1) * log(2 * pi) - 0.5 * log(Ti) - 0.5 * (Ti - 1) * log(sv2) -
      0.5 * as.numeric(t(et) %*% Si %*% et) + 0.5 * (ms^2 / ss2 - mu^2 / su2) +
      0.5 * log(ss2) + pnorm(ms / sqrt(ss2), log.p = TRUE) -
      (0.5 * log(su2) + pnorm(mu / sqrt(su2), log.p = TRUE))
  }
  tot
}

test_that("FD reports its exact log-likelihood at the optimum", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 5, N = 100, rand = 406, sig_u = 1, sig_v = 0.3,
    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))
  f <- suppressWarnings(psfm(y_fd ~ x_fd | z_fd, model_name = "FD", data = d,
    individual = "name", time = "year"))
  expect_equal(-f$opt$value, .fd_exact_loglik(f$out[, "par"], d), tolerance = 1e-8)
})
