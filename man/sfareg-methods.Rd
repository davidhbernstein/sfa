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
Standard extractor methods for objects of class \code{"sfareg"}, as
returned by \code{\link{psfm}}, \code{\link{sfm}}, \code{\link{zsfm}}, and
\code{\link{ttsfm}}. These let \code{"sfareg"} fits be used with the usual
R modeling idioms (e.g. \code{coef(fit)}, \code{logLik(fit)},
\code{AIC(fit)}, \code{BIC(fit)}) instead of reaching into the returned
list by hand.

\code{AIC()} and \code{BIC()} are not defined explicitly for
\code{"sfareg"} objects; they fall through to their default methods in the
\pkg{stats} package, which only require a working \code{logLik} method
(and, for \code{BIC()}, a working \code{nobs} method) to work correctly.
}
\usage{
\method{coef}{sfareg}(object, \dots)
\method{vcov}{sfareg}(object, type = c("hessian", "bhhh"), \dots)
\method{logLik}{sfareg}(object, \dots)
\method{nobs}{sfareg}(object, \dots)
}
\arguments{
\item{type}{For \code{vcov}: \code{"hessian"} (the default) for the inverse Hessian, or \code{"bhhh"} for the inverse outer product of the per-observation scores. See \sQuote{Details}.}

  \item{object}{An object of class \code{"sfareg"}.}
  \item{\dots}{Currently unused; included for S3 method consistency.}
}
\details{
\bold{BHHH standard errors.} \code{vcov(fit, type = "bhhh")} returns the inverse of the outer product of the per-observation scores instead of the inverse Hessian. Its value is that it needs no Hessian at all: it is defined whenever the scores are, which is exactly the case where the default path fails. When the Hessian is singular \code{vcov} otherwise falls back to a \emph{diagonal} approximation built from the reported standard errors, discarding every covariance; BHHH keeps them. It requires \code{keep_objective = TRUE} for the same reason \code{estfun} does, and is returned on the same scale \code{coef} reports (see below). It is not the default, because it is only as good as the information-matrix equality and so disagrees with the Hessian under misspecification -- which is useful when the two are compared on purpose and misleading when they are not.

\bold{Robust and clustered standard errors.} \code{vcov} returns the inverse Hessian, which is valid only if the likelihood is correctly specified and the observations are independent -- neither being safe to assume in applied frontier work. Registering \code{bread} and \code{estfun} makes the \pkg{sandwich} package operate on these fits, so \code{sandwich::sandwich()} gives Huber--White errors and \code{sandwich::vcovCL()} gives clustered ones, and either can be handed to \code{lmtest::coeftest()}:

\preformatted{fit <- sfm(y ~ x1 + x2, model_name = "NHN", data = d,
           keep_objective = TRUE)
lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = d$firm))}

Two requirements. The fit must have been made with \code{keep_objective = TRUE}, since the scores are differenced from the retained likelihood; and the \code{robust} divergence estimators are excluded, because they do not maximise a log-likelihood and so have no score. For \code{\link{psfm}} panel fits the scores are per \emph{firm}, so BHHH and sandwich standard errors are clustered by firm, and \code{bread()} scales by the number of firms rather than by \code{nobs()}. \code{sandwich::vcovHC()} does \emph{not} work and is not expected to: its finite-sample corrections are built from hat values, which are not defined for a nonlinear likelihood whose parameter vector includes variance components.

\bold{Which fits can supply scores.} \code{keep_objective = TRUE} is available on every entry point whose likelihood is separable into per-observation contributions: \code{\link{sfm}}, \code{\link{zsfm}}, \code{\link{lcsfm}}, \code{\link{ttsfm}} (\code{"TTNE"} and \code{"TTHN"}), \code{\link{copsfm}}, \code{\link{ivsfm}} (\code{"IVLIML"} and \code{"IVCF"}), \code{\link{selsfm}}, and \code{\link{psfm}} for \code{"GTRE"} (\code{estimator = "sml"}), \code{"TRE"}, \code{"GTRE_Z"}, \code{"TRE_Z"}, \code{"GTRE_FML"}, \code{"TFE_WMLE"}, \code{"FD"}, \code{"PL80"}, \code{"PL80_MVTN"}, \code{"BC92"}, \code{"K1990"} and \code{"K1990modified"}, whose contributions are per \emph{firm}. All but \code{sfm} and \code{psfm}'s four simulated-ML models are new in 1.2.1. The estimators that maximise no log-likelihood have no score and say so rather than returning an unusable fit: \code{ttsfm("TTNLS")} is nonlinear least squares, \code{ivsfm("C2SLS")} is 2SLS with a corrected intercept, \code{psfm}'s \code{"TFE"}, \code{"SSFE"}, \code{"SSRE"}, \code{"SSCRE"}, \code{"CSS"}, \code{"LS"}, \code{"KSS"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are moment-based, FE/LSDV or sequential, and \code{\link{npsfm}} is nonparametric.

\bold{The scale the scores live on.} \code{estfun} differences the retained likelihood at \code{fit$opt$par}, which is not always what \code{coef} reports: \code{ttsfm}, \code{copsfm} and \code{ivsfm} estimate \emph{log} sigmas and report sigmas, \code{ivsfm} reports \eqn{\rho = t/\sqrt{1 + t't}}, and \code{ivsfm("IVLIML")} additionally estimates reduced-form parameters it never reports. \code{vcov} therefore carries the fit's own Jacobian from the estimation scale onto the reported one -- the same delta method the reported standard errors already use -- and marginalises any parameters that are estimated but not reported by taking the corresponding block of the \emph{inverse}, not the inverse of the block. \code{estfun} itself returns the scores on the estimation scale, with that map attached as attributes, so its column count can exceed \code{length(coef(fit))}. From 1.2.1 this correction is applied on the Hessian path too: before it, \code{sqrt(diag(vcov(f)))} on a \code{copsfm} fit disagreed with the standard errors the same fit printed by a factor of \eqn{1/\hat\sigma_u}, and on an \code{ivsfm("IVLIML")} fit \code{vcov} failed outright.

The scores are taken by central differences with a per-parameter step, there being no analytic gradient for most of these models. Where one does exist (\code{NHN}) the two agree to about \eqn{3 \times 10^{-7}}, and they agree with \code{numDeriv} to the same order.

\code{coef.sfareg} returns the named vector of estimated parameters
(\code{object$coefficients}).

\code{vcov.sfareg} returns the estimated variance-covariance matrix. If the
fit was estimated with \code{optHessian = TRUE} and the Hessian is
invertible, this is \code{solve(object$opt$hessian)}. Otherwise it falls
back to a diagonal matrix built from the reported standard errors
(\code{object$std.errors}), with a warning; if neither is available it
returns a matrix of \code{NA}s with a warning.

\code{logLik.sfareg} returns the maximized log-likelihood as an object of
class \code{"logLik"}, with \code{"df"} and \code{"nobs"} attributes set so
that \code{AIC()}/\code{BIC()} work automatically. Every estimator in this
package minimizes the negative summed log-likelihood internally, so
\code{logLik = -object$opt$value}. \strong{Exception:} \code{psfm()}'s
\code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} models are estimated by a
moment-based procedure and \code{"SSFE"} by fixed-effects (LSDV)
regression, not by maximum likelihood; none of the three carries an
\code{$opt} component, so calling \code{logLik()} (or, transitively,
\code{AIC()}/\code{BIC()}) on those fits returns \code{NA} with a warning
rather than erroring.

\code{nobs.sfareg} returns the number of observations used in estimation,
recovered from \code{object$data} when the fitting function stored it, or
by re-evaluating the \code{data} argument of the original call otherwise.
Returns \code{NA} if neither is available.
}
\value{
\code{coef.sfareg} returns a named numeric vector.

\code{vcov.sfareg} returns a numeric matrix.

\code{logLik.sfareg} returns an object of class \code{"logLik"}.

\code{nobs.sfareg} returns a single integer (or \code{NA}).
}
\seealso{\code{\link{psfm}}, \code{\link{sfm}}, \code{\link{zsfm}}, \code{\link{ttsfm}}, \code{\link{print.sfareg}}, \code{\link{summary.sfareg}}}
