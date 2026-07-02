\name{psfm_bootstrap}
\alias{psfm_bootstrap}
\title{Parametric Bootstrap for \code{psfm()} Panel Stochastic Frontier Models}
\description{
Performs a parametric bootstrap for panel stochastic frontier models
(TRE / TRE_Z / GTRE / GTRE_Z) fit with \code{\link[sfa]{psfm}}.
Frontier, inefficiency (\code{u(z)}), and inefficiency (\code{h(z)})
parameter blocks are read directly from the fitted model's three-part
\code{\link[Formula]{Formula}} (\code{y ~ x | z_u | z_h}). For each of
\code{BOOT} replications, new disturbances are drawn, a simulated response is
constructed from the fitted parameters, the model is re-estimated in parallel,
and bootstrap standard errors / t-values are returned for every parameter in
\code{psfm_object$out}.
}
\usage{
psfm_bootstrap(psfm_object,
               numCores,
               BOOT,
               individual,
               h_type        = c("auto", "none", "scalar", "parametric"),
               maxit.psoptim = 1000,
               seed_offset   = 0,
               write_back    = TRUE,
               pkgs          = c("sfa", "Formula", "pbapply"),
               inefdec,
               rand.gtre     = NULL,
               rand.psoptim  = NULL,
               maxit.bobyqa  = 1,
               maxit.optim   = 1)
}
\arguments{
  \item{psfm_object}{A fitted model object of class \code{"sfareg"} returned by
    \code{\link[sfa]{psfm}}. Must contain \code{$out}, \code{$data},
    \code{$formula}, \code{$U}, \code{$H}, and \code{$model_name}.}
  \item{numCores}{Integer. Number of cores for the parallel cluster created
    internally via \code{\link[parallel]{makeCluster}}.}
  \item{BOOT}{Integer. Number of bootstrap replications (e.g. 199, 399, 999).}
  \item{individual}{Character. Name of the individual/panel-ID column in
    \code{psfm_object$data} (e.g. \code{"firm.id"}).}
  \item{h_type}{One of \code{"auto"}, \code{"none"}, \code{"scalar"}, or
    \code{"parametric"}, controlling how the individual-specific,
    time-invariant inefficiency term \code{h} is simulated:
    \describe{
      \item{\code{"none"}}{The model has no \code{h} term at all (e.g. plain
        TRE).}
      \item{\code{"scalar"}}{\code{h_i ~ abs(N(0, sigma_h))}, a single scalar
        standard deviation taken from the row of \code{$out} immediately
        following the \code{z} (\code{u}) block.}
      \item{\code{"parametric"}}{\code{h_i ~ abs(N(0, exp(data_h \%*\% h_vec)))},
        using the third right-hand-side part of the model formula as the
        \code{h}-covariates.}
      \item{\code{"auto"}}{(Default) Resolves to \code{"parametric"} if the
        model formula has three right-hand-side parts, otherwise
        \code{"none"}.}
    }
  }
  \item{maxit.psoptim}{Integer. Passed through to \code{psfm()} as
    \code{maxit.psoptim} on each bootstrap re-estimation. Default \code{1000}.}
  \item{seed_offset}{Integer added to the bootstrap index \code{b} when calling
    \code{set.seed()}, allowing reproducible-but-distinct seeds across multiple
    bootstrap runs. Default \code{0}.}
  \item{write_back}{Logical. If \code{TRUE} (default), the returned list
    includes \code{model}, a copy of \code{psfm_object} with the bootstrap
    standard errors and t-values written into \code{$out[,2]} and
    \code{$out[,3]}. If \code{FALSE}, only the raw bootstrap results are
    returned.}
  \item{pkgs}{Character vector of package names to load on each cluster
    worker via \code{\link[parallel]{clusterCall}}. Default
    \code{c("sfa", "Formula", "pbapply")}.}
  \item{inefdec}{Logical, with no default; must be supplied. Passed through to
    \code{psfm()} on each bootstrap re-estimation. Should match the
    \code{inefdec} setting used to originally fit \code{psfm_object}.}
  \item{rand.gtre}{Passed through to \code{psfm()} as \code{rand.gtre} on each
    bootstrap re-estimation. Default \code{NULL}.}
  \item{rand.psoptim}{Passed through to \code{psfm()} as \code{rand.psoptim} on
    each bootstrap re-estimation. Default \code{NULL}.}
  \item{maxit.bobyqa}{Integer. Passed through to \code{psfm()} as
    \code{maxit.bobyqa} on each bootstrap re-estimation. Default \code{1}.}
  \item{maxit.optim}{Integer. Passed through to \code{psfm()} as
    \code{maxit.optim} on each bootstrap re-estimation. Default \code{1}.}
}
\details{
The row layout of \code{psfm_object$out} is assumed to be:
\preformatted{
  row 1         : sigv
  row 2         : sigr
  next Kx rows  : frontier (x) parameters, in formula right-hand-side part 1 order
  next Kz rows  : u-block z parameters, in formula right-hand-side part 2 order
  next Kh rows  : h-block z parameters, in formula right-hand-side part 3 order
                  (only present when \code{h_type = "parametric"}; one extra
                  row for a scalar sigma_h when \code{h_type = "scalar"})
}
\code{psfm_bootstrap} validates that \code{nrow(psfm_object$out)} matches this
expected layout for the resolved \code{h_type} before running any bootstrap
replications, and stops with an informative error if it does not. It
similarly validates up front that \code{length(psfm_object$H)} equals the
number of unique individuals implied by the \code{individual} column, since
\code{$H} is one value per individual rather than one value per observation;
a mismatch here is also caught before any replications run, rather than after
the (potentially long-running) bootstrap completes.

For each replication \code{b}, new disturbances are drawn -- \code{v} (classic
symmetric noise), \code{u} (half-normal with covariate-driven sigma from the
\code{z} block), \code{r} (an individual random effect, constant within
individual and repeated over time), and \code{h} (per \code{h_type}, as
described above) -- and a simulated response
\code{y* = X \%*\% beta_hat + v + u + r + h} is constructed. The model is then
re-estimated on the simulated data via \code{sfa::psfm()}, executed in
parallel across replications using \code{\link[parallel]{parLapply}} (or
\code{\link[pbapply]{pblapply}} for a live progress bar, when the
\pkg{pbapply} package is installed). From each re-estimated model, both the
time-varying efficiency scores (\code{$U}) and the time-invariant,
individual-specific inefficiency scores (\code{$H}) are collected across
replications into \code{boot_eff} and \code{boot_eff_h} respectively.

Re-estimation failures (e.g. optimizer errors) are caught per replication,
recorded with their error message, and set to \code{NA} in the output
matrices rather than stopping the bootstrap; a single \code{\link{warning}} at
the end lists every failed replication index together with its error message.
Bootstrap standard errors are computed with \code{na.rm = TRUE}, so they
remain well-defined even when some replications failed.

Requires the \pkg{parallel} and \pkg{Formula} packages, and requires
\pkg{sfa} (and any other packages named in \code{pkgs}) to be installed and
loadable on each cluster worker.
}
\value{
A list with components:
  \item{boot_par}{A \code{BOOT} x (\code{nrow(psfm_object$out)} + 2) matrix of
    bootstrapped parameter draws, with the final two columns \code{"loglik"}
    and \code{"hours"}. Rows corresponding to failed replications are
    \code{NA}.}
  \item{boot_eff}{A \code{BOOT} x \code{nrow(psfm_object$data)} matrix of
    bootstrapped time-varying efficiency-score (\code{U}) draws. Rows
    corresponding to failed replications are \code{NA}.}
  \item{boot_eff_h}{A \code{BOOT} x \code{length(psfm_object$H)} matrix of
    bootstrapped time-invariant, individual-specific inefficiency-score
    (\code{H}) draws -- one column per unique individual (not one column
    per observation; \code{$H} does not vary over time within an
    individual). Columns are named by the individual ID values. Rows
    corresponding to failed replications are \code{NA}.}
  \item{se}{Named numeric vector of bootstrap standard errors for every
    parameter in \code{psfm_object$out}.}
  \item{tval}{Named numeric vector of bootstrap t-values
    (\code{estimate / boot SE}) for every parameter in \code{psfm_object$out}.}
  \item{failures}{Integer vector of replication indices \code{b} for which
    re-estimation failed, or \code{NULL} if none failed.}
  \item{model}{Only present when \code{write_back = TRUE}. A copy of
    \code{psfm_object} with \code{$out[,2]} and \code{$out[,3]} replaced by
    the bootstrap standard errors and t-values.}
}
\author{
David H. Bernstein, Chris Parmeter
}
\seealso{
\code{\link[sfa]{psfm}}, \code{\link[parallel]{parLapply}},
\code{\link[pbapply]{pblapply}}, \code{\link[Formula]{Formula}}
}
\examples{
\dontrun{
## GTRE_ZZ is a fitted "sfareg" object returned by sfa::psfm(), e.g.:
## GTRE_ZZ <- sfa::psfm(lnC_P3 ~ lnY + lnP1_P3 | Age | Class ,
##                      model_name = "GTRE_Z", data = data,
##                      individual = "MFI.ID")

boot_results <- psfm_bootstrap(
  psfm_object = GTRE_ZZ,
  numCores    = 4,
  BOOT        = 199,
  individual  = "MFI.ID",
  h_type      = "auto",
  inefdec     = FALSE
)

## Bootstrap standard errors / t-values for every model parameter
boot_results$se
boot_results$tval

## Bootstrap draws of time-varying (U) and time-invariant (H)
## efficiency/inefficiency scores, one row per replication
boot_results$boot_eff
boot_results$boot_eff_h

## Fitted object with bootstrap SEs/t-values written into $out
boot_results$model$out
}
}
\keyword{models}
\keyword{htest}
