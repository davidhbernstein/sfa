\name{cols_sfm}
\alias{cols_sfm}
\title{Corrected OLS for the normal-half-normal frontier}
\usage{
cols_sfm(formula, data, inefdec = TRUE)
}
\arguments{
\item{formula}{A two-sided formula, as for \code{\link{sfm}}.}
\item{data}{A data frame.}
\item{inefdec}{\code{TRUE} for a production frontier (inefficiency lowers the
response), \code{FALSE} for a cost frontier.}
}
\value{
A list with \code{coefficients} and \code{std.errors} for the frontier,
\code{gamma}, \code{sigmaSq} and their standard errors, the implied
\code{sigma_u}, \code{sigma_v} and \code{lambda}, the OLS \code{residuals},
\code{nobs}, and Coelli's \code{wald} statistic.
}
\description{
Estimates the normal-half-normal frontier by corrected OLS: invert the
second and third moments of the OLS residuals for \eqn{\sigma_u} and
\eqn{\sigma_v}, then shift the intercept by \eqn{E[u] = \sigma_u\sqrt{2/\pi}}
so it describes the frontier rather than the mean.
}
\details{
\bold{The standard errors are the point of this function.} Most applications
of COLS report the ordinary OLS standard errors for \eqn{\gamma} and
\eqn{\sigma^2}, which are not appropriate: those parameters are non-linear
functions of the residual moments \eqn{m_2} and \eqn{m_3}, not regression
coefficients. Coelli (1995, Appendix 1, equations A15-A16) gives the correct
delta-method variance,
\deqn{V(\hat\gamma) = g_2^2 V(m_2) + g_3^2 V(m_3) + 2 g_2 g_3 C(m_2, m_3),}
with \eqn{V(m_2) = \mu_4 - \mu_2^2}, \eqn{V(m_3) = \mu_6 - \mu_3^2 - 6\mu_2\mu_4
+ 9\mu_2^3} and \eqn{C(m_2,m_3) = \mu_5 - 4\mu_2\mu_3}, all divided by
\eqn{n}, and the higher moments re-expressed in terms of the estimates
themselves. That is what is computed here.

The frontier coefficients other than the intercept carry their OLS standard
errors, which are correct for them: COLS leaves the slopes untouched.
}
\references{
Coelli, T. (1995). Estimators and hypothesis tests for a stochastic frontier
function: A Monte Carlo analysis. \emph{Journal of Productivity Analysis}
\bold{6}, 247--268.

Olson, J. A., Schmidt, P. and Waldman, D. M. (1980). A Monte Carlo study of
estimators of stochastic frontier production functions. \emph{Journal of
Econometrics} \bold{13}, 67--82.
}
\seealso{\code{\link{sfm}}, \code{\link{inefficiency_test}}}
\examples{
d <- as.data.frame(data_gen_cs(N = 400, rand = 7, sig_u = 1, sig_v = 0.5,
                               cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                               a = 5, mu = 0.1))
cols_sfm(y_pcs ~ x1 + x2, data = d)
}
