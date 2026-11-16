\name{efficiency}
\alias{efficiency}
\title{Technical Efficiency Scores}
\description{
Extracts observation-specific technical efficiency from a fitted frontier,
choosing among the three usual point predictors and between the two scales the
dependent variable may be on.
}
\usage{
efficiency(object, type = c("bc", "jlms", "mode"), logDepVar = TRUE,
           newdata = NULL, component = c("transient", "persistent"))
}
\arguments{
\item{object}{An \code{"sfareg"} fit.}
\item{type}{Which predictor of \eqn{u_i} to use. \code{"bc"} (the default) is
  Battese--Coelli, \eqn{E[\exp(-u_i) \mid \varepsilon_i]}; \code{"jlms"} is
  Jondrow \emph{et al.}, \eqn{\exp(-E[u_i \mid \varepsilon_i])}; \code{"mode"}
  uses the mode of the posterior of \eqn{u_i}.}
\item{logDepVar}{\code{TRUE} (the default) if the dependent variable is on a
  log scale, so that efficiency is \eqn{\exp(-u_i)}; \code{FALSE} if it is on a
  level scale, giving \eqn{1 - u_i / f(x_i)} for a production frontier and
  \eqn{f(x_i) / (f(x_i) + u_i)} for a cost frontier (\code{inefdec = FALSE}).}
\item{newdata}{Data to rebuild the fitted frontier from, needed only when
  \code{logDepVar = FALSE}. Supply it when the data used for the fit can no
  longer be recovered from the stored call.}
\item{component}{Which efficiency, for the four-component panel models that
  report two. \code{"transient"} (the default) is \eqn{E[\exp(-u_{it}) \mid
  \cdot]}, one value per observation; \code{"persistent"} is
  \eqn{E[\exp(-h_i) \mid \cdot]}, one value per \emph{firm}. Only
  \code{psfm()}'s \code{"GTRE"}, \code{"GTRE_FML"} and \code{"GTRE_Z"} have
  a persistent component; asking any other fit for one is an error rather than
  a silent fall back to the transient score.}
}
\details{
\bold{Which predictor.} The three differ in what they report about the same
posterior, not in how well they are estimated.
\eqn{E[\exp(-u)] \ge \exp(-E[u])} by Jensen's inequality, so \code{"bc"} is
never below \code{"jlms"}; the gap widens as the posterior spreads out.
\code{"mode"} is the only one that can equal exactly 1, and it does so for every
observation whose posterior mean is negative -- typically much of the efficient
tail. That is a property of the predictor rather than a defect: where the
posterior mean of \eqn{u} is negative the single most likely value really is the
boundary. It is available only for the models whose posterior is a truncated
normal (\code{"NHN"}, \code{"NHN_Z"}, \code{"NE"}, \code{"NTN"}).

\bold{Which component.} A four-component panel fit has two efficiencies and
not one: a transient \eqn{u_{it}} that varies within firm and a persistent
\eqn{h_i} that does not. \code{efficiency()} returns the transient score by
default, matching \code{\link{predict.sfareg}}, and says so once per call
when the fit has both and the caller did not choose. \code{"TRE"} and
\code{"TRE_Z"} have a persistent random \emph{effect} but no persistent
\emph{inefficiency}, so they have no second score.

Note the lengths differ: a persistent score has one value per firm, so it
cannot be placed beside observation-level quantities without expanding it, and
\code{logDepVar = FALSE} is refused for it outright -- the fitted frontier is
observation-level and there is no ratio to form.

\bold{Which scale.} The package elsewhere assumes throughout that the
dependent variable is logged, which is the usual case and makes \eqn{u} a
proportional shortfall. \code{logDepVar = FALSE} instead forms
\eqn{(f(x_i) - u_i)/f(x_i)} for a production frontier, and for a cost frontier,
where \eqn{y = f(x) + u + v}, cost efficiency \eqn{f(x_i)/(f(x_i) + u_i)}:
minimum over actual cost, the level-scale counterpart of \eqn{\exp(-u_i)}. The
orientation is read from the fit's \code{inefdec}, including when it was passed
as a variable. That is a ratio of two estimated quantities and is
less well behaved than the log version: it is undefined where the fitted
frontier is zero and negative where \eqn{u} exceeds it, both of which happen
near zero output. A warning is issued if any score comes back negative or
non-finite, since silence there would be the real failure.
}
\value{
A numeric vector of efficiency scores: one per observation used in the fit, or
one per firm when \code{component = "persistent"}.
}
\seealso{
\code{\link{efficiency_ci}} for interval estimates of the same quantity.
}
\examples{
\donttest{
set.seed(2)
n  <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
y  <- 6 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.8))
d  <- data.frame(y = y, x1 = x1, x2 = x2)

fit <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)

summary(efficiency(fit))                      # Battese-Coelli, log scale
summary(efficiency(fit, type = "jlms"))       # never above the above
summary(efficiency(fit, type = "mode"))       # reaches exactly 1
summary(efficiency(fit, logDepVar = FALSE, newdata = d))
}
}
\keyword{models}
