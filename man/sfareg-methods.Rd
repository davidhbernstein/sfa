\name{sfareg-methods}
\alias{sfareg-methods}
\alias{coef.sfareg}
\alias{vcov.sfareg}
\alias{logLik.sfareg}
\alias{nobs.sfareg}
\alias{bread.sfareg}
\alias{estfun.sfareg}
\title{Extractor methods for "sfareg" objects}
\description{
Standard extractor methods for objects of class \code{"sfareg"}, which
\code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}},
\code{\link{ttsfm}}, \code{\link{lcsfm}}, \code{\link{copsfm}},
\code{\link{ivsfm}}, \code{\link{selsfm}} and \code{\link{esfm}} return.
They let a fit be used with the usual R idioms -- \code{coef(fit)},
\code{vcov(fit)}, \code{logLik(fit)}, \code{AIC(fit)}, \code{BIC(fit)} --
instead of reaching into the returned list by hand.

\tabular{ll}{
  \bold{method} \tab \bold{returns} \cr
  \code{coef} \tab the reported parameters, \code{object$coefficients} \cr
  \code{vcov} \tab their covariance: inverse Hessian (default), BHHH, sandwich or clustered \cr
  \code{logLik} \tab the maximized log-likelihood, with \code{"df"} and \code{"nobs"} for \code{AIC()} and \code{BIC()} \cr
  \code{nobs} \tab the number of observations the fit used \cr
  \code{estfun}, \code{bread} \tab scores and bread for the \pkg{sandwich} package (need \code{keep_objective = TRUE})
}

\code{AIC()} and \code{BIC()} are not defined for \code{"sfareg"}; the
\pkg{stats} defaults work through \code{logLik} and \code{nobs}.
}
\usage{
\method{coef}{sfareg}(object, \dots)
\method{vcov}{sfareg}(object, type = c("hessian", "bhhh", "sandwich", "clustered"),
  cluster = NULL, \dots)
\method{logLik}{sfareg}(object, \dots)
\method{nobs}{sfareg}(object, \dots)
}
\arguments{
\item{type}{For \code{vcov}: \code{"hessian"} (the default) for the inverse Hessian, \code{"bhhh"} for the inverse outer product of the per-observation scores, \code{"sandwich"} for the heteroskedasticity-robust sandwich, or \code{"clustered"} for its cluster-robust form. The last three need \code{keep_objective = TRUE}; the last two additionally need a Hessian, which \code{"bhhh"} does not. See \sQuote{Details}.}

\item{cluster}{For \code{vcov(type = "clustered")}: a vector or factor grouping the ROWS OF THE SCORE MATRIX. That is one row per observation for a cross-sectional fit and one row per FIRM for a panel one, so a firm-year vector is the wrong length by construction and is refused with a message saying so. A data frame is accepted for multiway clustering and is reduced with \code{interaction()}.}

  \item{object}{An object of class \code{"sfareg"}.}
  \item{\dots}{Currently unused; included for S3 method consistency.}
}
\details{
\subsection{What each method returns}{
\code{coef.sfareg} returns the named vector of reported parameters,
\code{object$coefficients}.

\code{vcov.sfareg} with the default \code{type = "hessian"} inverts the
Hessian stored on the fit (\code{object$opt$hessian}). Where the Hessian's
diagonal is positive it is first scaled to unit diagonal, so that parameters
on very different scales do not make an invertible matrix look singular.
Where the fit estimates on a different scale from the one it reports (see
\sQuote{The scale the scores live on} below), the inverse is then mapped onto
the reported scale. If the inversion or the mapping fails it warns and falls
back to a diagonal matrix built from the reported standard errors
(\code{object$std.errors}), which it also returns, without a warning, for a
fit that stores no Hessian. With neither it returns a matrix of \code{NA}s
with a warning. The other three types are described below.

\code{logLik.sfareg} returns the maximized log-likelihood as a
\code{"logLik"} object whose \code{"df"} and \code{"nobs"} attributes make
\code{AIC()} and \code{BIC()} work; see \sQuote{Two quantities that are easy
to confuse} for what \code{"df"} counts.

\code{nobs.sfareg} returns the number of observations the fit \emph{used},
not the rows supplied, so \code{BIC()} divides by the right \eqn{n}. It takes
the first of: the count stored on the fit (\code{object$nobs}); the rows of
stored data (\code{object$data}); the length of a per-observation vector the
fit stores (such as \code{exp_u_hat}); and, failing all of those, the rows of
the call's \code{data} argument re-evaluated in the caller's frames, accepted
only if it matches what the fit was made from. It returns \code{NA} if none of
these is available.
}

\subsection{Four covariances}{
\code{vcov(type = "sandwich")} returns the heteroskedasticity-robust sandwich and \code{vcov(type = "clustered")} its cluster-robust form. Both are formed on the ESTIMATION scale -- the scale the scores and the Hessian live on -- and only then mapped to the scale \code{coef} reports. That mapping is why they are not left to \code{sandwich::sandwich()}: that function composes bread and meat itself and has no way to know the two scales differ. For \code{ttsfm}, \code{copsfm}, \code{ivsfm} and \code{psfm}'s \code{PL80} family they do differ, and for \code{ivsfm("IVLIML")} they differ in dimension as well, since it maximises over 11 parameters and reports 6.

\code{sandwich::sandwich()} and \code{sandwich::vcovCL()} remain available and now return an estimation-scale covariance, which is internally consistent but not on the scale \code{coef} reports. Where a model reports what it estimates -- \code{sfm}, \code{zsfm}, \code{lcsfm}, \code{psfm}'s \code{TRE}/\code{GTRE} family -- the two routes agree exactly, and the cluster adjustment matches \code{vcovCL}'s default of \eqn{M/(M-1)}.
}

\subsection{BHHH standard errors}{
\code{vcov(fit, type = "bhhh")} returns the inverse of the outer product of the per-observation scores instead of the inverse Hessian. Its value is that it needs no Hessian at all: it is defined whenever the scores are, which is exactly the case where the default path fails. When the Hessian is singular \code{vcov} otherwise falls back to a \emph{diagonal} approximation built from the reported standard errors, discarding every covariance; BHHH keeps them. It requires \code{keep_objective = TRUE} for the same reason \code{estfun} does, and is returned on the same scale \code{coef} reports (see below). It is not the default, because it is only as good as the information-matrix equality and so disagrees with the Hessian under misspecification -- which is useful when the two are compared on purpose and misleading when they are not. At a parameter bound the scores in the constrained directions are not defined: the likelihood returns a stability penalty rather than a value, so \code{estfun} falls back to a one-sided difference on the admissible side, uses zero where neither side is admissible, and warns. BHHH errors built from such a fit understate the uncertainty in those directions.
}

\subsection{Robust and clustered standard errors with the sandwich package}{
\code{vcov} returns the inverse Hessian, which is valid only if the likelihood is correctly specified and the observations are independent -- neither being safe to assume in applied frontier work. Registering \code{bread} and \code{estfun} makes the \pkg{sandwich} package operate on these fits, so \code{sandwich::sandwich()} gives Huber--White errors and \code{sandwich::vcovCL()} gives clustered ones, and either can be handed to \code{lmtest::coeftest()}:

\preformatted{fit <- sfm(y ~ x1 + x2, model_name = "NHN", data = d,
           keep_objective = TRUE)
lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = d$firm))}

Two requirements. The fit must have been made with \code{keep_objective = TRUE}, since the scores are differenced from the retained likelihood; and the \code{robust} divergence estimators are excluded, because they do not maximise a log-likelihood and so have no score. For \code{\link{psfm}} panel fits the scores are per \emph{firm}, so BHHH and sandwich standard errors are clustered by firm, and \code{bread()} scales by the number of firms rather than by \code{nobs()}. \code{sandwich::vcovHC()} does \emph{not} work and is not expected to: its finite-sample corrections are built from hat values, which are not defined for a nonlinear likelihood whose parameter vector includes variance components.
}

\subsection{Which fits can supply scores}{
\code{keep_objective = TRUE} is available on every entry point whose likelihood is separable into per-observation contributions: \code{\link{sfm}}, \code{\link{zsfm}}, \code{\link{lcsfm}}, \code{\link{ttsfm}} (\code{"TTNE"} and \code{"TTHN"}), \code{\link{copsfm}}, \code{\link{ivsfm}} (\code{"IVLIML"} and \code{"IVCF"}), \code{\link{selsfm}}, and \code{\link{psfm}} for \code{"GTRE"} (\code{estimator = "sml"}), \code{"TRE"}, \code{"GTRE_Z"}, \code{"TRE_Z"}, \code{"GTRE_FML"}, \code{"TFE_WMLE"}, \code{"FD"}, \code{"PL80"}, \code{"PL80_MVTN"}, \code{"BC92"}, \code{"K1990"} and \code{"K1990modified"}, whose contributions are per \emph{firm}. All but \code{sfm} and \code{psfm}'s four simulated-ML models are new in 1.2.1. The estimators that maximise no log-likelihood have no score and say so rather than returning an unusable fit: \code{ttsfm("TTNLS")} is nonlinear least squares, \code{ivsfm("C2SLS")} is 2SLS with a corrected intercept, \code{psfm}'s \code{"TFE"}, \code{"SSFE"}, \code{"SSRE"}, \code{"SSCRE"}, \code{"CSS"}, \code{"LS"}, \code{"KSS"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are moment-based, FE/LSDV or sequential, and \code{\link{npsfm}} is nonparametric.
}

\subsection{The scale the scores live on}{
\code{estfun} differences the retained likelihood at \code{fit$opt$par}, which is not always what \code{coef} reports: \code{ttsfm}, \code{copsfm} and \code{ivsfm} estimate \emph{log} sigmas and report sigmas, \code{ivsfm} reports \eqn{\rho = t/\sqrt{1 + t't}}, and \code{ivsfm("IVLIML")} additionally estimates reduced-form parameters it never reports. \code{vcov} therefore carries the fit's own Jacobian from the estimation scale onto the reported one -- the same delta method the reported standard errors already use -- and marginalises any parameters that are estimated but not reported by taking the corresponding block of the \emph{inverse}, not the inverse of the block. \code{estfun} itself returns the scores on the estimation scale, with that map attached as attributes, so its column count can exceed \code{length(coef(fit))}. From 1.2.1 this correction is applied on the Hessian path too: before it, \code{sqrt(diag(vcov(f)))} on a \code{copsfm} fit disagreed with the standard errors the same fit printed by a factor of \eqn{1/\hat\sigma_u}, and on an \code{ivsfm("IVLIML")} fit \code{vcov} failed outright.

The scores are taken by central differences with a per-parameter step, there being no analytic gradient for most of these models. Where one does exist (\code{NHN}) the two agree to about \eqn{3 \times 10^{-7}}, and they agree with \code{numDeriv} to the same order.
}

\subsection{Two quantities that are easy to confuse}{
\code{object$opt$value} is
the OBJECTIVE the optimizer minimized, which is the negative summed
log-likelihood for most but not all estimators here, and
\code{object$coefficients} is what the fit REPORTS, which is not always what
it estimated. \code{logLik()} therefore reads neither blindly: \code{"df"} is
the number of parameters ESTIMATED (\code{length(object$opt$par)} where that
exists), because \code{\link{ivsfm}}\code{("IVLIML")} maximizes over 11 and
reports 6; and the returned value is the PLAIN log-likelihood, because
\code{\link{lcsfm}} with \code{penalty_c > 0} maximizes a penalized
likelihood and carries the unpenalized one separately.
}

\subsection{Estimators with no log-likelihood}{
\code{\link{sfm}}'s robust
divergence estimators (\code{robust = "mlqe"}, \code{"psi"} or
\code{"mdpd"}) minimize a divergence rather than a negative log-likelihood,
so \code{logLik()} returns \code{NA} with a warning and \code{AIC()}/
\code{BIC()} propagate it. Their objective is not on the likelihood scale and
must not be compared with one; \code{print()} names it accordingly.
\strong{Exception:} \code{psfm()}'s
\code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} models are estimated by a
moment-based procedure and \code{"SSFE"} by fixed-effects (LSDV)
regression, not by maximum likelihood; none of the three carries an
\code{$opt} component, so calling \code{logLik()} (or, transitively,
\code{AIC()}/\code{BIC()}) on those fits returns \code{NA} with a warning
rather than erroring.
}
}
\value{
\code{coef.sfareg} returns a named numeric vector.

\code{vcov.sfareg} returns a numeric matrix.

\code{logLik.sfareg} returns an object of class \code{"logLik"}.

\code{nobs.sfareg} returns a single integer (or \code{NA}).
}
\seealso{\code{\link{print.sfareg}}, \code{\link{summary.sfareg}}, \code{\link{TIC}}, \code{\link{influence_sfa}}; the entry points listed under \sQuote{Description}.}
