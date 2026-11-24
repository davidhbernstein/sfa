## Marginal effects of the variance determinants on inefficiency

## phi(a)/Phi(a), the inverse Mills ratio, computed in logs so the far left
## tail does not become 0/0. Used by the truncated-normal marginal effect.
.sfa_lambda_ratio <- function(a) {
  exp(stats::dnorm(a, log = TRUE) - stats::pnorm(a, log.p = TRUE))
}

## The pre-truncation mean, per observation, rebuilt from the block the fit
## stores rather than read off a `$mu` component. There is no such component:
## `object$mu` PARTIAL-MATCHES `mu_spec` and silently hands back the spec list,
## which is the same trap the A48 note records for `object$data`. Everything
## below therefore uses [[ ]] exact indexing.
##
## het.R builds mu as `Zmu %*% delta` -- an IDENTITY link, NOT the `link` field
## that mu_spec carries beside it (that field describes the SCALE's link and
## must not be applied here). A homoskedastic truncated normal keeps its single
## mu in `out` instead.
.sfa_me_mu <- function(object, n) {
  ms <- object[["mu_spec"]]
  m <- NULL
  if (!is.null(ms) && !is.null(ms[["Z"]]) && !is.null(ms[["delta"]]) &&
    ncol(ms[["Z"]]) == length(ms[["delta"]])) {
    m <- as.numeric(ms[["Z"]] %*% ms[["delta"]])
  }
  if (is.null(m)) {
    o <- object[["out"]]
    if (!is.null(o) && "mu" %in% rownames(o)) m <- o["mu", "par"]
  }
  if (is.null(m) || !length(m) || !all(is.finite(m))) {
    stop("marginal_effects(): this truncated-normal fit carries no usable ",
      "pre-truncation mean, so dE[u]/dz cannot be formed. A fit made with ",
      "muhet = ~ z stores one; see ?marginal_effects.",
      call. = FALSE
    )
  }
  rep_len(m, n)
}

## delta values for `nm`, in that order, 0 where the design has no such column.
.sfa_me_delta_by_name <- function(spec, nm) {
  out <- stats::setNames(rep(0, length(nm)), nm)
  if (is.null(spec) || is.null(spec[["delta"]])) {
    return(out)
  }
  hit <- intersect(nm, names(spec[["delta"]]))
  out[hit] <- as.numeric(spec[["delta"]][hit])
  out
}



## Columns of a design that are 0/1 dummies (exactly the values 0 and 1).
.sfa_me_is_binary <- function(x) {
  u <- unique(x)
  length(u) == 2L && all(sort(u) == c(0, 1))
}

## Groups of dummy columns that model.matrix() made from ONE factor, read off
## the design's "assign" and "contrasts" attributes. Only designs built from a
## formula argument (uhet, muhet, scaling) keep them; a pipe segment's design
## does not, and its dummies are then switched one at a time. A group must be
## all 0/1 and mutually exclusive, so an interaction with a factor is not
## mistaken for one. Returns a list of column-name vectors of length >= 2.
.sfa_me_factor_groups <- function(M) {
  asg <- attr(M, "assign")
  ctr <- attr(M, "contrasts")
  if (is.null(M) || is.null(asg) || is.null(ctr) || length(asg) != ncol(M)) {
    return(list())
  }
  cn <- colnames(M)
  out <- list()
  for (fn in names(ctr)) {
    for (a in unique(asg[startsWith(cn, fn)])) {
      cols <- cn[asg == a]
      if (length(cols) < 2L || !all(startsWith(cols, fn))) next
      sub <- M[, cols, drop = FALSE]
      if (!all(apply(sub, 2, .sfa_me_is_binary)) || any(rowSums(sub) > 1)) next
      out[[length(out) + 1L]] <- cols
    }
  }
  out
}

## E[u] and Var[u] per observation with some columns of the designs SET to
## given values, everything else at its observed value. `set` is a named
## numeric vector; a name is set in every design that has that column (the
## scale, the pre-truncation mean, the scaling factor), so a covariate in
## several designs moves through all of them. This is the discrete-change
## counterpart of the derivatives below: the moments are recomputed from the
## full formula at each value, not by a scale shortcut, so the truncated
## normal's mu and sigma_u move together.
.sfa_me_moments_at <- function(zs, ss, ms, mu_fixed, scaling, set = NULL) {
  put <- function(M) {
    if (is.null(M) || is.null(set)) return(M)
    for (nm in intersect(names(set), colnames(M))) M[, nm] <- set[[nm]]
    M
  }
  eta <- as.numeric(put(zs[["Z"]]) %*% zs[["delta"]])
  sig <- if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))
  h <- 1
  if (scaling) {
    h <- exp(as.numeric(put(ss[["Z"]]) %*% ss[["delta"]]))
    sig <- sig * h
  }
  if (identical(zs$family, "halfnormal")) {
    return(list(e = sig * sqrt(2 / pi), v = sig^2 * (1 - 2 / pi)))
  }
  if (identical(zs$family, "exponential")) {
    return(list(e = sig, v = sig^2))
  }
  mu <- if (!is.null(ms) && !is.null(ms[["Z"]]) && !is.null(ms[["delta"]]) &&
    ncol(ms[["Z"]]) == length(ms[["delta"]])) {
    as.numeric(put(ms[["Z"]]) %*% ms[["delta"]])
  } else {
    mu_fixed
  }
  mu <- rep_len(mu, length(sig)) * h
  a <- mu / sig
  lam <- .sfa_lambda_ratio(a)
  list(e = mu + sig * lam, v = sig^2 * (1 - lam * (lam + a)))
}

## Constant columns of the z design are dropped: the derivative with respect
## to an intercept is not a marginal effect. Discrete determinants get a
## discrete change, not a derivative (`discrete`; see ?marginal_effects).
marginal_effects <- function(object, average = FALSE, component = c("u", "h"),
                             discrete = "auto") {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit, as returned by sfm() or psfm().",
      call. = FALSE
    )
  }
  ## "h" selects the PERSISTENT block of psfm("GTRE_Z"); every other model has
  ## only the one.
  component <- match.arg(component)
  if (identical(component, "h")) {
    if (is.null(object$z_spec_h)) {
      stop(
        "component = \"h\" is only available for a psfm(model_name = \"GTRE_Z\") ",
        "fit with a third formula segment, y ~ x | z | zp, which is what ",
        "parameterizes the persistent sigma_h. This fit is model_name \"",
        if (is.null(object$model_name)) "unknown" else object$model_name,
        "\".",
        call. = FALSE
      )
    }
    zs <- object$z_spec_h
  } else {
    zs <- object$z_spec
  }
  if (is.null(zs)) {
    stop(
      "marginal_effects() needs a model whose inefficiency scale depends on ",
      "covariates. This fit is model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\", which has a single homoskedastic sigma_u and therefore no ",
      "marginal effect to report. Use one of the _Z models with a formula ",
      "such as y ~ x1 + x2 | z.",
      call. = FALSE
    )
  }

  Z <- zs[["Z"]]
  delta <- zs[["delta"]]
  eta <- as.numeric(Z %*% delta)

  ## sigma_u on the scale each family actually uses.
  sigma_u <- if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))

  ## Under the scaling property (sfm(scaling = ~ z), truncated normal only) the
  ## covariates drive ONE factor h = exp(z'delta_s) that multiplies both the
  ## scale and the pre-truncation mean; Zu and Zmu are intercept-only there, so
  ## the sigma_u just computed is the scalar sigma_u0 and misses h entirely.
  ## Take h on board, and differentiate with respect to the SCALING design.
  ss <- object[["s_spec"]]
  scaling <- identical(component, "u") && !is.null(ss) &&
    !is.null(ss[["Z"]]) && ncol(ss[["Z"]]) > 0L
  h <- 1
  if (scaling) {
    h <- exp(as.numeric(ss[["Z"]] %*% ss[["delta"]]))
    sigma_u <- sigma_u * h
    Z <- ss[["Z"]]
    delta <- ss[["delta"]]
  }

  ## Moments of u given sigma_u.
  if (identical(zs$family, "halfnormal")) {
    e_u <- sigma_u * sqrt(2 / pi)
    var_u <- sigma_u^2 * (1 - 2 / pi)
  } else if (identical(zs$family, "exponential")) {
    e_u <- sigma_u
    var_u <- sigma_u^2
  } else if (identical(zs$family, "truncnormal")) {
    ## mu enters E[u] separately from the scale, so E[u] is NOT proportional to
    ## sigma_u and the outer(e_u, half * d) form below does not apply (gap M4).
    ## The pre-truncation mean comes off the fit rather than being rebuilt.
    mu_i <- .sfa_me_mu(object, length(sigma_u)) * h
    a <- mu_i / sigma_u
    lam <- .sfa_lambda_ratio(a)
    g <- 1 - lam * (lam + a)
    e_u <- mu_i + sigma_u * lam
    var_u <- sigma_u^2 * g
  } else {
    stop("unsupported inefficiency family: ", zs$family, call. = FALSE)
  }

  ## d sigma_u / d z_k = scale_k * sigma_u, with scale_k = delta_k under the
  ## SD link and delta_k/2 under the variance link.
  half <- if (identical(zs$link, "sd")) 1 else 0.5

  varying <- function(M) {
    if (is.null(M) || !ncol(M)) return(character(0))
    colnames(M)[vapply(seq_len(ncol(M)), function(j) length(unique(M[, j])) > 1L,
                       logical(1))]
  }
  nm <- varying(Z)
  ## The truncated normal has a second design, for mu. A covariate may sit in
  ## it alone (muhet = ~ z2 beside uhet = ~ z1, or muhet with no uhet), and it
  ## still moves E[u]; taking the names from Z alone used to drop it.
  if (!scaling && identical(zs$family, "truncnormal")) {
    nm <- union(nm, varying(object[["mu_spec"]][["Z"]]))
  }
  if (!length(nm)) {
    stop("every column of the variance-determinant design is constant, so ",
      "there is no covariate to differentiate with respect to.",
      call. = FALSE
    )
  }
  ## delta for each reported covariate, 0 where the scale design lacks it.
  d <- .sfa_me_delta_by_name(list(delta = stats::setNames(delta, colnames(Z))), nm)

  if (scaling) {
    ## a = mu/sigma_u is CONSTANT in z here, because h cancels out of it, so
    ## E[u] = h E[u*] and Var[u] = h^2 Var[u*] whatever the family is. The
    ## effect is therefore exact and needs no family branch at all, and
    ## delta_s_k is itself the semi-elasticity d log E[u] / d z_k.
    me_e <- outer(e_u, d)
    me_v <- outer(var_u, 2 * d)
  } else if (identical(zs$family, "truncnormal")) {
    ## Chain rule through BOTH designs. d sigma_u/d z_k = half * delta_u_k *
    ## sigma_u as for the other families; d mu/d z_k = delta_mu_k, because mu
    ## carries an IDENTITY link (het.R builds it as Zmu %*% delta), NOT the
    ## z_link that mu_spec records beside it. A covariate may sit in one design,
    ## the other, or both, so the two delta vectors are matched BY NAME and a
    ## covariate absent from a design contributes zero through it.
    dmu <- .sfa_me_delta_by_name(object[["mu_spec"]], nm)
    ## lambda'(a) = -lambda (lambda + a), so
    ##   dE/dmu    = 1 - lambda (lambda + a) = g
    ##   dE/dsigma = lambda + a lambda (lambda + a)
    ## and for Var = sigma^2 g(a), with g'(a) = lambda (lambda + a)(2 lambda + a) - lambda.
    dE_dmu <- g
    dE_dsig <- lam + a * lam * (lam + a)
    gp <- lam * (lam + a) * (2 * lam + a) - lam
    me_e <- matrix(0, length(sigma_u), length(nm), dimnames = list(NULL, nm))
    me_v <- me_e
    for (k in seq_along(nm)) {
      ds <- half * d[[k]] * sigma_u
      dm <- dmu[[k]]
      da <- dm / sigma_u - a * half * d[[k]]
      me_e[, k] <- dE_dmu * dm + dE_dsig * ds
      me_v[, k] <- 2 * sigma_u * ds * g + sigma_u^2 * gp * da
    }
  } else {
    me_e <- outer(e_u, half * d)
    me_v <- outer(var_u, 2 * half * d)
  }
  ## Name the columns after the component actually differentiated, so a u table
  ## and an h table cannot be confused once separated from their call.
  colnames(me_e) <- paste0("dE_", component, ".d", nm)
  colnames(me_v) <- paste0("dVar_", component, ".d", nm)
  type <- stats::setNames(rep("derivative", length(nm)), nm)

  ## ---- discrete changes --------------------------------------------------
  ## A dummy has no derivative to speak of: E[u] at z_k = 1 minus E[u] at
  ## z_k = 0, per observation, other covariates at their observed values. A
  ## named non-binary column (an ordinal coded 1, 2, 3, entering the index
  ## linearly) gets the change per unit across its observed range,
  ## [E(max) - E(min)] / (max - min), so every observation is used. Dummies
  ## from one factor() are switched as a group: each level against the
  ## reference, its sibling dummies held at 0.
  disc <- .sfa_me_resolve_discrete(discrete, nm, object, zs, ss, scaling)
  if (length(disc$cols)) {
    ## A homoskedastic mu has no design to set; it is the fitted constant.
    ms <- if (identical(zs$family, "truncnormal")) object[["mu_spec"]] else NULL
    mu_fixed <- if (identical(zs$family, "truncnormal") && (is.null(ms) || is.null(ms[["Z"]]))) {
      .sfa_me_mu(object, length(sigma_u))
    } else {
      NA_real_
    }
    mom <- function(set) .sfa_me_moments_at(zs, ss, ms, mu_fixed, scaling, set)
    for (k in disc$cols) {
      j <- match(k, nm)
      grp <- disc$group[[k]]
      if (!is.null(grp)) {
        on <- stats::setNames(as.numeric(grp == k), grp)
        off <- stats::setNames(rep(0, length(grp)), grp)
        hi <- mom(on)
        lo <- mom(off)
        step <- 1
        type[[k]] <- "factor level"
      } else {
        x <- disc$values[[k]]
        hi <- mom(stats::setNames(max(x), k))
        lo <- mom(stats::setNames(min(x), k))
        step <- max(x) - min(x)
        type[[k]] <- if (.sfa_me_is_binary(x)) "discrete" else "one step"
      }
      me_e[, j] <- (hi$e - lo$e) / step
      me_v[, j] <- (hi$v - lo$v) / step
      colnames(me_e)[j] <- paste0("DeltaE_", component, ".", k)
      colnames(me_v)[j] <- paste0("DeltaVar_", component, ".", k)
    }
  }
  names(type) <- NULL

  out <- data.frame(sigma_u, e_u, var_u, me_e, me_v, check.names = FALSE)
  names(out)[1:3] <- paste0(c("sigma_", "E_", "Var_"), component)

  avg <- colMeans(cbind(me_e, me_v))
  attr(out, "average") <- avg
  attr(out, "link") <- zs$link
  attr(out, "family") <- zs$family
  attr(out, "component") <- component
  attr(out, "effect_type") <- stats::setNames(rep(type, 2L), colnames(cbind(me_e, me_v)))

  if (isTRUE(average)) {
    return(avg)
  }
  out
}


## Which reported covariates get a discrete change, and how.
##   discrete = NULL / FALSE   none: every column is a derivative (pre-1.3.0)
##   discrete = "auto"         every 0/1 column, plus factor groups
##   discrete = c("a", "b")    exactly these columns (0/1 or ordinal);
##                             include "auto" to add the 0/1 columns too
## Returns cols (in nm order), the observed values of each (from whichever
## design holds it), and the factor group a column belongs to, if any.
.sfa_me_resolve_discrete <- function(discrete, nm, object, zs, ss, scaling) {
  none <- list(cols = character(0), values = list(), group = list())
  if (is.null(discrete) || isFALSE(discrete)) return(none)
  if (!is.character(discrete) || !length(discrete) || anyNA(discrete)) {
    stop("`discrete` must be \"auto\", a character vector of determinant ",
      "names, or NULL for derivatives throughout.",
      call. = FALSE
    )
  }
  designs <- if (scaling) {
    list(ss[["Z"]])
  } else {
    list(zs[["Z"]], if (identical(zs$family, "truncnormal")) object[["mu_spec"]][["Z"]])
  }
  designs <- Filter(Negate(is.null), designs)
  values <- list()
  for (k in nm) {
    for (M in designs) {
      if (k %in% colnames(M)) {
        values[[k]] <- M[, k]
        break
      }
    }
  }
  named <- setdiff(discrete, "auto")
  bad <- setdiff(named, nm)
  if (length(bad)) {
    stop("`discrete` names ", paste(dQuote(bad, FALSE), collapse = ", "),
      ", which ", if (length(bad) == 1L) "is" else "are", " not among this ",
      "fit's non-constant determinants: ", paste(dQuote(nm, FALSE), collapse = ", "), ".",
      call. = FALSE
    )
  }
  auto <- if ("auto" %in% discrete) {
    nm[vapply(nm, function(k) .sfa_me_is_binary(values[[k]]), logical(1))]
  } else {
    character(0)
  }
  cols <- nm[nm %in% union(auto, named)]
  group <- list()
  for (M in designs) {
    for (g in .sfa_me_factor_groups(M)) {
      for (k in intersect(g, cols)) group[[k]] <- g
    }
  }
  list(cols = cols, values = values, group = group)
}
