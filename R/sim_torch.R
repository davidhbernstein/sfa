## Optional torch evaluator for the GTRE and TRE simulated likelihoods.
## It computes the same quantity as .gtre_sim_logdens() over the same
## firm-specific draws, on a CPU or GPU through the torch package (Suggests).
## The optimizer and every default are unchanged; see
## notes/code_history/sim_torch.md.

## Resolve backend, device and precision once per fit.
.sim_torch_setup <- function(sim_backend, sim_device, sim_dtype, model_name) {
  if (identical(sim_backend, "R")) {
    if (!is.null(sim_device) || !is.null(sim_dtype)) {
      warning("`sim_device` and `sim_dtype` apply only with sim_backend = ",
        "\"torch\"; ignored.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (!model_name %in% c("GTRE", "TRE")) {
    stop("sim_backend = \"torch\" is implemented for model_name \"TRE\" and ",
      "for \"GTRE\" with estimator = \"sml\" only, not ", dQuote(model_name),
      ". Those are the two simulated likelihoods with a stacked (all firms ",
      "and draws at once) implementation for torch to evaluate.",
      call. = FALSE
    )
  }
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("sim_backend = \"torch\" needs the 'torch' package: ",
      "install.packages(\"torch\"), then torch::install_torch().",
      call. = FALSE
    )
  }
  if (!isTRUE(torch::torch_is_installed())) {
    stop("The 'torch' package is installed but its libtorch backend is not. ",
      "Run torch::install_torch() once, then restart R.",
      call. = FALSE
    )
  }

  ## Device: the fastest available unless one is named.
  if (is.null(sim_device)) {
    sim_device <- if (isTRUE(torch::cuda_is_available())) {
      "cuda"
    } else if (isTRUE(torch::backends_mps_is_available())) {
      "mps"
    } else {
      "cpu"
    }
  }
  if (!is.character(sim_device) || length(sim_device) != 1L ||
    !sim_device %in% c("cpu", "cuda", "mps")) {
    stop("`sim_device` must be \"cpu\", \"cuda\" or \"mps\".", call. = FALSE)
  }
  if (identical(sim_device, "cuda") && !isTRUE(torch::cuda_is_available())) {
    stop("sim_device = \"cuda\" but torch finds no CUDA device.", call. = FALSE)
  }
  if (identical(sim_device, "mps") && !isTRUE(torch::backends_mps_is_available())) {
    stop("sim_device = \"mps\" but torch finds no Apple GPU (MPS) device.",
      call. = FALSE
    )
  }

  ## Precision: double wherever the device has it. Apple's MPS has no float64,
  ## so there the only choice is single precision, and the caller is told.
  if (is.null(sim_dtype)) {
    sim_dtype <- if (identical(sim_device, "mps")) "float32" else "float64"
  }
  if (!is.character(sim_dtype) || length(sim_dtype) != 1L ||
    !sim_dtype %in% c("float64", "float32")) {
    stop("`sim_dtype` must be \"float64\" or \"float32\".", call. = FALSE)
  }
  if (identical(sim_device, "mps") && identical(sim_dtype, "float64")) {
    stop("Apple's MPS device has no float64. Use sim_dtype = \"float32\", or ",
      "sim_device = \"cpu\" for double precision.",
      call. = FALSE
    )
  }
  list(
    device = sim_device, dtype_name = sim_dtype,
    dtype = if (identical(sim_dtype, "float64")) torch::torch_float64() else torch::torch_float32()
  )
}

## Move the fixed pieces of the stacked likelihood to the device once.
## gid is the firm index of each stacked row (1-based, as torch's R interface
## uses for index_add_).
.sim_torch_prep <- function(setup, yv, Xall, H1, H2, gid, N) {
  dev <- setup$device
  dt <- setup$dtype
  tt <- function(a) torch::torch_tensor(a, dtype = dt, device = dev)
  list(
    setup = setup,
    yv = tt(as.numeric(yv)),
    X = tt(unname(as.matrix(Xall))),
    H1 = tt(unname(H1)),
    H2 = if (is.null(H2)) NULL else tt(unname(H2)),
    gid = torch::torch_tensor(as.integer(gid), dtype = torch::torch_long(), device = dev),
    N = as.integer(N),
    R = ncol(H1),
    chunks = .sim_torch_chunks(length(yv), ncol(H1)),
    tt = tt
  )
}

## Row blocks of at most `cells` observation-by-draw cells. The whole n x R
## array at once fell off a cliff on a 16 GB Apple GPU: at 1e7 cells one
## evaluation took 26 s against 4 s for the R path, while at 2.5e6 it was 8x
## faster. Blocks accumulate into the per-firm sums by index_add, so any row
## split gives the same likelihood. Same idea as SFrontiers.jl's `chunks`.
.sim_torch_chunks <- function(n, R, cells = getOption("sfa.torch_chunk_cells", 1e6)) {
  rows <- max(1L, as.integer(floor(cells / R)))
  starts <- seq.int(1L, n, by = rows)
  lapply(starts, function(s) c(s, min(n, s + rows - 1L)))
}

## log Phi(z), stable in single precision. In double, log(erfc(-z/sqrt 2)/2)
## is exact down to the -37 clip. In single precision Phi underflows below
## about z = -13, so below -10 this switches to the Mills-ratio series
## log Phi(z) = -z^2/2 - log(-z) - log(2 pi)/2 + log(1 - z^-2 + 3 z^-4
## - 15 z^-6 + 105 z^-8), whose truncation error at z = -10 is about 1e-8,
## below single precision.
.torch_log_pnorm <- function(z, single) {
  direct <- torch::torch_log(torch::torch_erfc(-z / sqrt(2)) / 2)
  if (!single) return(direct)
  zs <- torch::torch_clamp(z, max = -10)
  w <- torch::torch_reciprocal(zs * zs)
  series <- (zs * zs) * (-0.5) - torch::torch_log(-zs) - 0.5 * log(2 * pi) +
    torch::torch_log(w * (-1) + 1 + (w^2) * 3 - (w^3) * 15 + (w^4) * 105)
  torch::torch_where(z < -10, series, direct)
}

## Per-firm simulated log-likelihood, the torch counterpart of
## .gtre_sim_logdens(): E is the n x R matrix of composed errors (u - v
## orientation already applied). Returns an N-vector of doubles in R.
.gtre_sim_logdens_torch <- function(E, lambda, sigma, prep) {
  lt <- .torch_log_terms(E, lambda, sigma, identical(prep$setup$dtype_name, "float32"))
  lp <- torch::torch_zeros(prep$N, prep$R, dtype = prep$setup$dtype, device = prep$setup$device)
  lp <- lp$index_add(1L, prep$gid, lt)
  torch::torch_logsumexp(lp, dim = 2L) - log(prep$R)
}

## log of (2/sigma) phi(E/sigma) Phi(-lambda E/sigma), elementwise.
.torch_log_terms <- function(E, lambda, sigma, single) {
  C <- .SFA_CONSTANTS
  z1 <- torch::torch_clamp(E / sigma, min = C$CLIP_Z1_LOWER, max = C$CLIP_Z1_UPPER)
  z2 <- torch::torch_clamp(-E * (lambda / sigma), min = C$CLIP_Z2_LOWER, max = C$CLIP_Z2_UPPER)
  ## Every expression puts the TENSOR first. In R's torch, `scalar - tensor`
  ## converts the scalar to a float32 tensor before the arithmetic, so a
  ## leading R double silently costs a float64 evaluation about 1e-8 per
  ## term. `tensor - scalar` keeps the tensor's precision.
  (z1 * z1) * (-0.5) + .torch_log_pnorm(z2, single) +
    (log(2 / sigma) - 0.5 * log(2 * pi))
}

## The whole stacked likelihood for one parameter vector. Same layout,
## orientation and mixture as fn_1_loop() in psfm(): x = (lambda, sigma,
## sigma_r[, sigma_h], beta). Returns per-firm log-likelihoods (doubles).
.gtre_sim_ll_torch <- function(x, prep, model_name, inefdec_n, n_x_vars) {
  gtre <- identical(model_name, "GTRE")
  k0 <- if (gtre) 4L else 3L
  beta <- prep$tt(x[(k0 + 1L):(k0 + n_x_vars)])
  lam <- x[1]
  sig <- x[2]
  single <- identical(prep$setup$dtype_name, "float32")
  zeros <- function() torch::torch_zeros(prep$N, prep$R, dtype = prep$setup$dtype,
                                         device = prep$setup$device)
  ## inefdec_n * (y - X b): the production-oriented residual.
  base_all <- (prep$yv - prep$X$matmul(beta)) * inefdec_n
  lp1 <- zeros()
  lp2 <- if (gtre) zeros() else NULL
  for (ch in prep$chunks) {
    s <- ch[1]
    e <- ch[2]
    b <- base_all[s:e]$unsqueeze(2L)
    g <- prep$gid[s:e]
    if (gtre) {
      ## alpha enters as -/+ inefdec_n * sigma_r * draw (both halves of the
      ## symmetric mixture), h as + sigma_h * |draw|.
      a <- prep$H1[s:e, ] * (inefdec_n * x[3])
      h <- prep$H2[s:e, ] * x[4]
      lp1 <- lp1$index_add(1L, g, .torch_log_terms(b - a + h, lam, sig, single))
      lp2 <- lp2$index_add(1L, g, .torch_log_terms(b + a + h, lam, sig, single))
    } else {
      ## TRE: alpha follows the orientation (A108), so the composed error is
      ## inefdec_n * (y - X b) - sigma_r * draw.
      lp1 <- lp1$index_add(1L, g, .torch_log_terms(b - prep$H1[s:e, ] * x[3], lam, sig, single))
    }
  }
  l1 <- torch::torch_logsumexp(lp1, dim = 2L) - log(prep$R)
  ll <- if (gtre) {
    l2 <- torch::torch_logsumexp(lp2, dim = 2L) - log(prep$R)
    torch::torch_logsumexp(torch::torch_stack(list(l1, l2), dim = 2L), dim = 2L) - log(2)
  } else {
    l1
  }
  as.numeric(ll$cpu()$to(dtype = torch::torch_float64()))
}
