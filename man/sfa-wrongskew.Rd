\name{sfa-wrongskew}
\alias{sfa-wrongskew}
\alias{wrongskew}
\title{The wrong-skewness suite: what to do when the residuals lean the wrong way}
\description{
A guide to the options this package offers when the OLS residuals of a production frontier are skewed positively -- the \dQuote{wrong} way -- and the estimated inefficiency collapses to zero. This page is documentation only; each option is a function or an argument documented in its own right.
}
\details{
\bold{The problem.} Waldman (1982) showed that for the Normal-Half Normal frontier, if the third central moment of the OLS residuals is positive then the OLS fit is a local maximum of the likelihood at \eqn{\hat\sigma_u = 0}. The fit then reports an industry in which every firm is fully efficient. Horrace and Wright (2020) extend the result to the truncated-normal and doubly-truncated-normal cases. The same thing happens to the moment estimator for a different reason: \eqn{m_3 = \sigma_u^3\sqrt{2/\pi}(1-4/\pi)} has no real solution with \eqn{\sigma_u > 0} when \eqn{m_3 > 0}. Olson, Schmidt and Waldman (1980) call these the Type I failure (no admissible \eqn{\sigma_u}) and the Type II failure (no admissible \eqn{\sigma_v}).

\bold{It is not a data problem.} Simar and Wilson (2009) make the point that has to be kept in view: observing wrong skew is entirely consistent with the model being correctly specified. It happens routinely when \eqn{\sigma_u} is small relative to \eqn{\sigma_v}, and the smaller the sample the more often. So respecifying the model in response to it is, in their phrase, a solution in search of a problem -- unless there is an independent reason to prefer the new specification.

\bold{What the package offers.} They divide into three kinds.

\emph{Keep the model, change the estimator.}
\describe{
  \item{\code{\link{sfm}(estimator = "acols")}}{Parmeter and Zhao (2023). Corrected OLS with the first absolute moment \eqn{E|\varepsilon|} in place of the third central moment. The third moment is the only part of the inversion that can fail, so removing it removes the failure: on the authors' own design at \eqn{\lambda = 0.4}, \eqn{n = 200}, Type I failures fall from 48\% to 15\%.}
  \item{\code{\link{sfm}(estimator = "cmle")}}{Zhao and Parmeter (2022). Maximum likelihood subject to two moment conditions the assumed density implies. The Waldman stationary point does not satisfy them, so the search is kept away from it. On wrong-skew samples the proportion of fits returning \eqn{\hat\sigma_u \approx 0} falls from about 96\% to about 3\%.}
}

\emph{Keep the model, change the question.}
\describe{
  \item{\code{\link{skewness_test}}, \code{\link{symmetry_test}}}{Ask first whether the residual skew is significantly wrong at all, rather than reading its sign off a point estimate. Often it is not, and nothing further is needed.}
  \item{\code{\link{skewness_decomp}}}{Bonanno, De Giovanni and Domma (2017). Splits the composed third moment into the parts due to the inefficiency, to the noise, and to the dependence between them, and reports Zenga's median-based measure alongside. Under the textbook specification two of the three are identically zero, which is the whole reason a positive residual skew has nowhere to go.}
  \item{\code{\link{moment_range}}}{What skewness a given pair of distributions can produce at all. A residual skew outside the achievable range is evidence about the distributions, not about the firms.}
}

\emph{Change the model.}
\describe{
  \item{\code{\link{esfm}}}{Hafner, Manner and Simar (2018). Widens the distribution of \eqn{u} by a finite-sample correction that admits either sign of skew and converges to the assumed density as \eqn{n} grows.}
  \item{\code{\link{copsfm}}}{Bonanno, De Giovanni and Domma (2017), generalizing Smith (2008). A generalized-logistic noise, which can itself be skewed either way, and a copula linking \eqn{u} and \eqn{v}. Either ingredient can account for a positive composed skew without setting \eqn{\sigma_u} to zero.}
  \item{\code{\link{sfm}(model_name = "NB")}}{Carree (2002). A binomial inefficiency, the one distribution here that can be skewed either way: \eqn{p > 1/2} means most firms carry considerable inefficiency and only a few sit near the frontier. Same sample moment, opposite economics. Corrected OLS only.}
  \item{\code{\link{sfm}(model_name = "NGB2")}}{A generalized beta of the second kind inefficiency, whose shape parameters accommodate a much wider range of third and fourth moments than the half normal.}
}

\bold{Choosing.} If the skew is not significant, stop. If it is, prefer an estimator change to a model change unless there is a substantive reason for the new distribution, since a respecified model answers a different question. \code{"cmle"} and \code{"acols"} are the cheapest options and leave the model alone. Among the respecifications, \code{\link{copsfm}} is the most informative when the mechanism matters, because \code{\link{skewness_decomp}} then says which ingredient is doing the work.
}
\references{
Bonanno, G., De Giovanni, D. and Domma, F. (2017). The \dQuote{wrong skewness} problem: a re-specification of stochastic frontiers. \emph{Journal of Productivity Analysis} 47(1), 49-64. \doi{10.1007/s11123-017-0492-8}

Carree, M. A. (2002). Technological inefficiency and the skewness of the error component in stochastic frontier analysis. \emph{Economics Letters} 77(1), 101-107. \doi{10.1016/S0165-1765(02)00119-2}

Hafner, C. M., Manner, H. and Simar, L. (2018). The \dQuote{wrong skewness} problem in stochastic frontier models: A new approach. \emph{Econometric Reviews} 37(4), 380-400. \doi{10.1080/07474938.2016.1140284}

Olson, J. A., Schmidt, P. and Waldman, D. M. (1980). A Monte Carlo study of estimators of stochastic frontier production functions. \emph{Journal of Econometrics} 13(1), 67-82. \doi{10.1016/0304-4076(80)90043-3}

Parmeter, C. F. and Zhao, S. (2023). An alternative corrected ordinary least squares estimator for the stochastic frontier model. \emph{Empirical Economics} 64(6), 2831-2857. \doi{10.1007/s00181-023-02401-1}

Simar, L. and Wilson, P. W. (2009). Inferences from cross-sectional, stochastic frontier models. \emph{Econometric Reviews} 29(1), 62-98. \doi{10.1080/07474930903324523}

Waldman, D. M. (1982). A stationary point for the stochastic frontier likelihood. \emph{Journal of Econometrics} 18(2), 275-279. \doi{10.1016/0304-4076(82)90041-0}

Zhao, S. and Parmeter, C. F. (2022). The \dQuote{wrong skewness} problem: Moment constrained maximum likelihood estimation of the stochastic frontier model. \emph{Economics Letters} 221, 110901. \doi{10.1016/j.econlet.2022.110901}
}
\seealso{\code{\link{sfm}}, \code{\link{copsfm}}, \code{\link{esfm}}, \code{\link{skewness_decomp}}, \code{\link{skewness_test}}, \code{\link{moment_range}}}
\keyword{documentation}
