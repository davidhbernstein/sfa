\name{tt_boundary_report}
\alias{tt_boundary_report}
\alias{print.sfa_tt_boundary}
\title{Is a collapsed two-tier noise scale the maximum, or did the optimizer stop short?}
\usage{
tt_boundary_report(object, data, sigma_v_starts = c(0.05, 0.3, 1),
                   tol = 1e-3, tol_loglik = 1e-4, ...)
}
\arguments{
\item{object}{A \code{\link{ttsfm}} fit with \code{model_name = "TTHN"} or \code{"TTNE"}.}
\item{data}{The data the model was fitted to. \code{ttsfm()} does not retain it and every refit here needs it.}
\item{sigma_v_starts}{Positive starting values for \eqn{\sigma_v}. One refit is run from each, with every other parameter started at the fitted values.}
\item{tol}{Below this, the returned \eqn{\sigma_v} counts as collapsed.}
\item{tol_loglik}{How much higher a restart's log-likelihood has to be before it counts as higher at all. Restarts landing on the same optimum differ in the last few digits; without this they would be reported as an improvement of zero log units.}
\item{...}{Passed to \code{\link{ttsfm}} for each refit.}
}
\description{
Refits a two-tier frontier from several starting values of \eqn{\sigma_v}, scores every solution on the same likelihood, and reports which is highest. Use it when \code{\link{ttsfm}} returns \eqn{\hat\sigma_v \approx 0} and you need to know whether that is the maximum likelihood estimate or a failure of the optimizer.
}
\details{
\bold{This is a likelihood comparison and not a test.} No p-value is reported and none should be inferred. The hypothesis a reader would naturally attach to it -- \eqn{\sigma_v = 0} -- sits on the boundary of the parameter space, where the usual chi-square reference distribution for a likelihood-ratio statistic does not apply, and no appropriate null distribution has been established for this model. The output is a table of log-likelihoods and a statement about which is largest.

\bold{Why it is needed.} \code{"TTHN"} drives \eqn{\sigma_v} to zero when \eqn{\sigma_u} and \eqn{\sigma_w} are close, and that collapse is the maximum likelihood estimate rather than a defect: with \eqn{\sigma_u = \sigma_w} the difference \eqn{w - u} is already symmetric and heavier-tailed than a normal, so a small normal \eqn{v} underneath it is very nearly invisible. On a 2,000-observation draw at \eqn{\sigma_u = \sigma_w = 1} the profile in \eqn{\sigma_v} is flat to 0.000 log units from \eqn{10^{-6}} upward and decreasing above 0.01, with no interior maximum. Distinguishing that from an optimizer that merely stopped early is exactly what this function does.

The ordinary normal linear model on the same data is reported alongside, because a collapsed \eqn{\sigma_v} is not the same claim as \dQuote{there is no frontier here}: on the design above the collapsed two-tier fit still beats the linear model by 9.6 log units, with both tiers estimated at 1.081 against truths of 1.0.

\bold{Cost.} One \code{ttsfm()} refit per entry in \code{sigma_v_starts}. That is inexpensive for \code{"TTNE"} and expensive for \code{"TTHN"}, where a single fit at 2,000 observations can take several minutes.
}
\value{
An object of class \code{"sfa_tt_boundary"}, a list with components:
\item{table}{A data frame with one row per start: the start, the \eqn{\sigma_v} it reached, and the log-likelihood.}
\item{best}{Which row attained the highest log-likelihood. Ties go to the returned fit.}
\item{returned_is_best}{\code{TRUE} when no restart beat the fit as returned by more than \code{tol_loglik}.}
\item{gap}{The best restart's log-likelihood less the returned fit's.}
\item{collapsed}{\code{TRUE} when the returned \eqn{\sigma_v} is below \code{tol}.}
\item{lm_logLik, lm_gap}{The ordinary linear model's log-likelihood, and how far the fit beats it.}
}
\references{
Polachek, S. and Yoon, B. J. (1987). A two-tiered earnings frontier estimation of employer and employee information in the labor market. \emph{The Review of Economics and Statistics} 69(2), 296-302. \doi{10.2307/1927239}
}
\seealso{\code{\link{ttsfm}}, and \code{\link{sfa_diagnostics}} for optimizer diagnostics generally.}
\examples{
\donttest{
d <- data_gen_cs(N = 1200, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                 beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0, sig_w = 0.5)
f <- ttsfm(y_ttne ~ x1 + x2, data = d, model_name = "TTNE")
tt_boundary_report(f, data = d, sigma_v_starts = c(0.05, 0.8))
}
}
\keyword{models}
