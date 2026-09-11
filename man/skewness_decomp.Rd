\name{skewness_decomp}
\alias{skewness_decomp}
\alias{print.sfa_skew_decomp}
\title{Where the skewness of the composed error comes from}
\usage{
skewness_decomp(object, nodes = 160L)
}
\arguments{
\item{object}{A fitted model of class \code{"sfareg"}: a \code{\link{copsfm}} fit, or an \code{\link{sfm}} fit with \code{model_name = "NHN"} or \code{"NE"}. Other models error, because the decomposition needs the two error components to have separately identified third moments.}

\item{nodes}{Number of Gauss-Legendre nodes per dimension for the cross moments and the median. The default is where the FGM case stops moving in the ninth significant digit; raise it if a copula with strong tail dependence is in use.}
}
\description{
Splits the third central moment of the composed error \eqn{\varepsilon = v - u} into the part contributed by the asymmetry of the inefficiency term, the part contributed by the asymmetry of the noise, and the part contributed by the dependence between them, following Bonanno, De Giovanni and Domma (2017).
}
\details{
Write \eqn{\tilde u = u - E[u]} and \eqn{\tilde v = v - E[v]}. Then
\deqn{E[(\varepsilon - E\varepsilon)^3] = E[\tilde v^3] - E[\tilde u^3] - 3E[\tilde v^2 \tilde u] + 3E[\tilde v \tilde u^2],}
and the last two terms vanish under independence. The function reports the three pieces as \code{inefficiency} (\eqn{-E[\tilde u^3]}), \code{noise} (\eqn{E[\tilde v^3]}) and \code{dependence} (the two cross terms together).

The point of the decomposition is what it says about the textbook model. A normal \eqn{v} sets the noise term to zero, and independence sets the dependence term to zero, so under the standard specification the sign of the composed third moment \emph{is} the sign of \eqn{-E[\tilde u^3]} and nothing else. Every one-sided inefficiency distribution in common use is positively skewed, so a production frontier forces a negative composed skew, and a positively skewed residual has nowhere to go but \eqn{\hat\sigma_u = 0}. Skew in the noise and dependence between the components are perfectly ordinary explanations for the same residual, and the textbook model cannot express either. That is the sense in which Bonanno, De Giovanni and Domma call the wrong-skewness anomaly an ill-posed problem.

Also reported is Zenga's median-based measure, their Eq. (13),
\deqn{E[(\varepsilon - Me(\varepsilon))^3] = E[(\varepsilon - E\varepsilon)^3] + [E\varepsilon - Me(\varepsilon)]^3 + 3[E\varepsilon - Me(\varepsilon)]\mathrm{Var}(\varepsilon),}
whose sign is the one that says which way the density actually leans. When it disagrees in sign with the third central moment, \code{print()} says so: the mean-based measure is the misleading one in that case.

The marginal contributions are closed forms. The dependence term is computed by Gauss-Legendre quadrature on the unit square against the copula density, which covers every marginal and copula \code{copsfm()} supports rather than only the exponential/generalized-logistic/FGM combination for which the paper gives an explicit formula. The two agree to the quadrature's own accuracy. Note that the closed forms printed as Eq. (10) and Eq. (11) of the paper do not: their dependence terms disagree with a direct evaluation of the same expectations, and the corrected expressions -- with the Gauss hypergeometric weights \eqn{A_1 = -\delta_v[\psi(2\alpha_v) - \psi(\alpha_v)]} and \eqn{A_2 = \delta_v^2[\psi'(\alpha_v) - \psi'(2\alpha_v) - (\psi(2\alpha_v)-\psi(\alpha_v))^2]} -- are the ones pinned in \code{tests/testthat/test-wrongskew-suite.R}.
}
\value{
An object of class \code{"sfa_skew_decomp"}, a list with components:
\item{components}{Named numeric of length three: \code{inefficiency}, \code{noise}, \code{dependence}. They sum to \code{total}.}
\item{total}{The third central moment of the composed error implied by the fit.}
\item{share}{\code{components / total}.}
\item{mean, median, variance}{The first moment, median and variance of the composed error implied by the fit.}
\item{median_based}{Zenga's measure, Eq. (13).}
\item{model_name, udist, vdist, copula, copula_par}{What the decomposition was taken of.}
}
\references{
Bonanno, G., De Giovanni, D. and Domma, F. (2017). The \dQuote{wrong skewness} problem: a re-specification of stochastic frontiers. \emph{Journal of Productivity Analysis} 47(1), 49-64. \doi{10.1007/s11123-017-0492-8}

Zenga, M. (1985). \emph{Statistica descrittiva}. Turin: Giappichelli.
}
\seealso{\link{sfa-wrongskew} for the other ways this package handles wrong skewness, \code{\link{copsfm}} for the model the decomposition was written for, \code{\link{skewness_test}} for a formal test on the residuals.}
\examples{
set.seed(21)
n <- 300
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- 0.5 + 0.8 * d$x1 - 0.4 * d$x2 + rnorm(n, 0, 0.5) - abs(rnorm(n, 0, 1))

## Under the textbook model two of the three components are identically zero.
f <- sfm(y ~ x1 + x2, data = d, model_name = "NHN", estimator = "cols")
skewness_decomp(f)

\donttest{
## With a skewed noise and a copula, all three are live.
g <- copsfm(y ~ x1 + x2, data = d, copula = "fgm",
            udist = "exponential", vdist = "glogistic", n_nodes = 64)
skewness_decomp(g)
}
}
\keyword{models}
