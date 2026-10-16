\name{sfa-troubleshooting}
\alias{sfa-troubleshooting}
\alias{troubleshooting}
\title{When a fit warns: what the message means and what to do about it}
\description{
One page for the warnings, \code{NOTE}s and \code{NA}s the fitting functions produce, with the symptom first and the remedy after it. This page is documentation only; everything it points at is a function or an argument documented in its own right. For the specific case of positively skewed residuals and a zero \eqn{\sigma_u}, see \code{\link{sfa-wrongskew}}, which covers that literature on its own.
}
\details{
\bold{Run this first.} \code{\link{sfa_diagnostics}(fit)} reports how the optimizer actually finished -- the convergence code and message, the eigenvalues and conditioning of the Hessian, the implied parameter correlations, and, if the fit was made with \code{keep_objective = TRUE}, the gradient at the reported optimum and a likelihood slice through each parameter. Most of the questions below are answered by its output. \code{print} and \code{summary} of a fit also print a \code{NOTE} for each boundary state the fit reached, so read those before reading the coefficients.

\bold{A boundary, or a failure? They are different, and the package distinguishes them.} A one-sided scale of zero is sometimes the maximum likelihood estimate and sometimes a sign that the fit went wrong, and which one it is can be read off the residual skew. The flags \code{$wrong_skew}, \code{$sigma_u_at_bound} and \code{$residual_m3} exist to be read together.

\describe{
  \item{\dQuote{\code{sigma_u} is on the zero boundary and the OLS residuals have the wrong skew}}{The boundary \emph{is} the maximum likelihood estimate: with a non-negative third central moment the likelihood has no interior maximum (the Type I failure of Olson, Schmidt and Waldman 1980; Waldman 1982). This is not a numerical failure and refitting will not move it. Read it as no evidence of inefficiency in these data, treat the efficiency scores as uninformative, and use \code{\link{skewness_test}} to ask whether the skew is significantly wrong at all before concluding anything. \code{\link{sfa-wrongskew}} lists what the package offers if it is.}
  \item{\dQuote{\code{sigma_u} is on the zero boundary, but the OLS residuals are skewed the right way}}{The opposite case, and the suspect one. The data say inefficiency is present, so the boundary is unlikely to be the maximum and the fit may not be one. Refit from other starting values with \code{start_val}, and inspect \code{$opt$convergence}, \code{$wrong_skew}, \code{$sigma_u_at_bound} and \code{$residual_m3}. The Waldman reading above does \emph{not} apply here.}
  \item{\code{sfm(model_name = "NTN")}: \dQuote{\code{sigma_v} has collapsed to the boundary}}{\code{lambda} has run away, so the fit attributes all the variation to inefficiency and none to noise. On such a sample the likelihood has no interior maximum in \code{lambda}: it rises monotonically toward a finite limit as \code{lambda} grows. The boundary is therefore the estimate rather than a failure, the likelihood is flat there, and \emph{every} standard error is \code{NA} because no finite one exists at that point. Read it as \dQuote{no noise component identified in these data}, and see \code{$ntn_sigma_v_at_bound}.}
  \item{\code{sfm(model_name = "tHN")}: \code{sigma_u} collapsed}{Heavy-tailed noise absorbing the one-sided component is a known property of this model rather than a numerical failure, so this one has its own flag, \code{$thn_sigma_u_at_bound}, and its own reading. Treat the efficiency scores as uninformative when it is set.}
  \item{\code{tHN}: \dQuote{inspect \code{$thn_starts}}}{\code{tHN} is fitted from several starting values. \code{$thn_starts$n_distinct > 1} means the likelihood is multimodal for these data and the returned fit is only the best optimum \emph{found}. Profiling over a grid of \code{nu} is the honest next step; \code{nu} is weakly identified near the normal limit, where a large \code{nu} and a near-singular Hessian are expected rather than surprising.}
}

\bold{The optimizer.}
\describe{
  \item{\dQuote{The final \code{optim()} stage did not produce a usable result}}{\code{$opt$convergence} is \code{99}, which is this package's own code: the last of the three minimizer stages returned a non-finite objective or a non-finite Hessian, so the estimates come from the preceding stage and the standard errors are \code{NA}. The estimates themselves are usable. Try \code{optHessian = FALSE}, \code{\link[=vcov.sfareg]{vcov}(type = "bhhh")}, or different \code{start_val}.}
  \item{\dQuote{non-finite value supplied by \code{optim}}}{The finite-difference step left the region where the likelihood is representable. This is most common when a scale is already close to zero. Different starting values usually avoid it; so does fitting the nested simpler model first and starting from its estimates.}
  \item{A convergence code of 1, 51 or 52}{\code{1} is \code{optim}'s iteration limit -- raise \code{maxit.optim}, and \code{maxit.bobyqa}, \code{maxit.nlminb} or \code{maxit.psoptim} for the earlier stages. \code{51} and \code{52} are \code{L-BFGS-B}'s own warning and error reports from its line search, with the detail in \code{$opt$message}. Neither is automatically fatal and neither is automatically harmless, so \code{\link{sfa_diagnostics}} does not guess: it calls such a code \emph{benign} only when it can also show a small relative gradient and a positive definite Hessian at the reported point, and \emph{unverified} when the fit carries no objective to check them with. Refit with \code{keep_objective = TRUE} to get the verdict rather than the guess.}
}

\bold{Standard errors.}
\describe{
  \item{All of them are \code{NA}}{The Hessian is singular or non-finite at the reported point. If the fit is at one of the boundaries above, that is the answer and no standard error exists there. Otherwise \code{\link[=vcov.sfareg]{vcov}(type = "bhhh")} needs no Hessian and is often defined where the Hessian is not; it inverts the outer product of the scores instead, and says so explicitly when that is singular too. \code{\link{sfa_diagnostics}} reports which parameters the likelihood cannot identify.}
  \item{\code{vcov(type = "sandwich")} or \code{"clustered"} refuses}{Both need the per-observation scores. Refit with \code{keep_objective = TRUE}; the error says so, and reports whether \code{"bhhh"} would be defined for the same fit.}
  \item{The simulated-ML models}{\code{sfm} with \code{"NLN"} or \code{"NW"}, and \code{\link{psfm}}'s \code{"GTRE"} family, maximise a \emph{simulated} likelihood, so part of the reported uncertainty is simulation noise. \code{\link{simulation_se}} measures how much by refitting the same data with independent draw sequences. If estimates move between runs, that is the quantity to report.}
}

\bold{Reading the output.}
\describe{
  \item{\code{logLik()}, \code{AIC()} or \code{BIC()} returns \code{NA} with a warning}{Several \code{\link{psfm}} estimators -- \code{"GTRE_SEQ1"}, \code{"GTRE_SEQ2"}, \code{"SSFE"} and the other moment-based and FE/LSDV ones -- maximise no likelihood at all, so no log-likelihood exists to return. Use them for the estimates, not for a likelihood-based comparison.}
  \item{A \code{delta} is twice, or half, what you expected}{The variance-determinant link differs by model family and this is a real difference rather than a typo. \code{sfm}'s \code{_Z} models put \eqn{z'\delta} on the standard deviation, \eqn{\sigma_u = \exp(z'\delta)}; \code{psfm}'s \code{TRE_Z} and \code{GTRE_Z} put it on the variance, as the competing packages do. Fitting the same data both ways gives the same \eqn{\sigma_u} and the same log-likelihood, but \eqn{\delta} itself differs by a factor of two. \code{sfm}'s \code{z_link = "var"} aligns the two.}
  \item{Every efficiency score is near 1, or they all look alike}{Check the boundary flags first: a collapsed \eqn{\sigma_u} makes the scores uninformative by construction, and no amount of reading them will say otherwise. Where the fit is sound, \code{\link{efficiency_ci}} gives Horrace and Schmidt (1996) intervals for each \eqn{u_i}, which is the right way to see how little or how much a point score is pinning down.}
  \item{The row names are not the parameters you expected}{Several models report a reparameterization rather than raw scales -- \code{lambda} and \code{sigma} instead of \code{sigv} and \code{sigu} for \code{NHN} and \code{NTN}, \code{theta} rather than \code{sigu} for \code{NU}. The returned \code{$out} is one row per parameter with columns \code{par}, \code{st_err} and \code{t-val}; index it as \code{fit$out[, "par"]}.}
}

\bold{What to change, in order.} Starting values (\code{start_val}) first, since most of the states above are about where the optimizer stopped. Then the iteration limits (\code{maxit.bobyqa}, \code{maxit.psoptim}, \code{maxit.optim}). Then the estimator, where the model offers one -- \code{\link{sfm}}'s \code{estimator = "acols"} or \code{"cmle"} exist precisely for the wrong-skew case. Only then the model itself: a different \code{model_name} is a different assumption about the data, not a numerical remedy, and respecifying in response to a boundary that is the correct answer is a solution in search of a problem (Simar and Wilson 2009).
}
\references{
Horrace, W. C. and Schmidt, P. (1996). Confidence statements for efficiency estimates from stochastic frontier models. \emph{Journal of Productivity Analysis} 7(2-3), 257-282. \doi{10.1007/BF00157044}

Olson, J. A., Schmidt, P. and Waldman, D. M. (1980). A Monte Carlo study of estimators of stochastic frontier production functions. \emph{Journal of Econometrics} 13(1), 67-82. \doi{10.1016/0304-4076(80)90043-3}

Simar, L. and Wilson, P. W. (2009). Inferences from cross-sectional, stochastic frontier models. \emph{Econometric Reviews} 29(1), 62-98. \doi{10.1080/07474930903324523}

Waldman, D. M. (1982). A stationary point for the stochastic frontier likelihood. \emph{Journal of Econometrics} 18(2), 275-279. \doi{10.1016/0304-4076(82)90041-0}
}
\seealso{\code{\link{sfa-wrongskew}}, \code{\link{sfa_diagnostics}}, \code{\link{skewness_test}}, \code{\link{efficiency_ci}}, \code{\link{simulation_se}}, \code{\link{sfareg-methods}}, \code{\link{sfm}}, \code{\link{psfm}}}
\keyword{documentation}
