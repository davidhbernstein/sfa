## Density-power weights (gap L17): they respond to an observation far from the fitted
## surface, not to a bad regressor the surface bends toward. See
## notes/code_history/robust_weights.md.

density_weights <- function(object, sigma_v = NULL, sigma_u = NULL, c = NULL,
                            normalize = TRUE) {
  if (is.numeric(object)) {
    e <- object
    if (is.null(sigma_v) || is.null(sigma_u) || is.null(c))
      stop("When 'object' is a residual vector, sigma_v, sigma_u and c are required.",
           call. = FALSE)
  } else {
    e       <- .robust_residuals(object)
    sigma_v <- sigma_v %||% .robust_get(object, "sigma_v")
    sigma_u <- sigma_u %||% .robust_get(object, "sigma_u")
    c       <- c       %||% .robust_get(object, "c")
    if (is.null(c)) stop("No tuning parameter found on the fit; supply 'c'.",
                         call. = FALSE)
  }
  if (!is.finite(c) || c < 0) stop("c must be non-negative.", call. = FALSE)
  if (c <= 1e-10) return(rep(1, length(e)))

  f <- .dens_nhn(e, sigma_v, sigma_u)
  w <- f^c
  if (isTRUE(normalize)) w <- w / max(w)
  w
}
