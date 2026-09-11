\name{endogeneity_test}
\alias{endogeneity_test}
\alias{print.sfa_endog_test}
\title{Were the regressors endogenous after all?}
\description{
A Wald test of \eqn{H_0\!: \rho = 0} on an \code{\link{ivsfm}} fit -- the hypothesis that the noise \eqn{v} is uncorrelated with the reduced-form errors of the endogenous regressors, so that the endogeneity correction was not needed and a plain \code{\link{sfm}} fit would have been consistent. This is the Durbin--Wu--Hausman idea carried into the stochastic frontier, and is Equation (25) of Hou, Ramalho and Roseta-Palma (2025).
}
\usage{endogeneity_test(object, level = 0.05)}
\arguments{
\item{object}{A fit from \code{\link{ivsfm}} with \code{model_name = "IVLIML"} or \code{"IVCF"}. \code{"C2SLS"} is a moment estimator and carries no correction parameter, so it is refused.}
\item{level}{Size used for the reported verdict.}
}
\value{
An object of class \code{c("sfa_endog_test", "htest")} carrying \code{statistic} (the Wald \eqn{W}), \code{parameter} (its degrees of freedom, the number of endogenous regressors), \code{p.value}, the estimated \code{rho} with its standard errors and z-ratios, and \code{reject}.
}
\details{
\bold{The statistic.} With \eqn{\hat\rho} the vector of correlations between \eqn{v} and the reduced-form errors,

\deqn{W = \hat\rho' \,[\widehat{\mathrm{Var}}(\hat\rho)]^{-1}\, \hat\rho \;\sim\; \chi^2(q),}

\eqn{q} the number of endogenous regressors. With a single endogenous regressor it is exactly the squared z-ratio already printed for \code{rho_} in the fit.

\bold{Why this is an ordinary chi-square.} \eqn{\rho} lives in the open unit ball, so \eqn{\rho = 0} is an \emph{interior} point of the parameter space. That is the opposite of the situation \code{\link{inefficiency_test}} faces, where \eqn{\sigma_u = 0} sits on a boundary and the null distribution is a chi-bar-square mixture. No mixture is needed here.

\bold{Standard errors for \eqn{\rho}.} \code{ivsfm} parameterises the correlations as \eqn{\rho = t/\sqrt{1 + t't}}, which keeps them inside the unit ball for any real \eqn{t}. Before version 1.2.0 the standard errors of \eqn{\rho} were reported as \code{NA} on the ground that this map has a non-diagonal Jacobian. It does, but the Jacobian is short --

\deqn{\partial \rho_i/\partial t_j = \delta_{ij}/s - t_i t_j / s^3, \qquad s = \sqrt{1 + t't},}

that is \eqn{J = (I - \rho\rho')/s} -- so \eqn{\mathrm{Var}(\hat\rho) = J\,\mathrm{Var}(\hat t)\,J'} by the delta method. At the null \eqn{J = I}, which is what makes the test well calibrated there. \code{ivsfm} now reports these standard errors and carries the full \eqn{\mathrm{Var}(\hat\rho)} as \code{$vcov_rho}, which is what the \emph{joint} test needs.

\bold{Size and power.} Measured in this package over 1,000 replications per cell on the design of Hou, Ramalho and Roseta-Palma's section 4.1.1, with \eqn{T = 20} and \eqn{N} varying. The \eqn{\rho = 0} block is size; the rest is power. Monte Carlo standard errors at 1,000 replications are 0.0095, 0.0069 and 0.0031 at the three levels.

\tabular{lrrrrr}{
  \eqn{\rho} \tab \eqn{n} \tab 10\% \tab 5\% \tab 1\% \tab \eqn{se/sd} \cr
  0.0 \tab 500 \tab 0.101 \tab 0.051 \tab 0.014 \tab 0.999 \cr
  0.0 \tab 1000 \tab 0.098 \tab 0.051 \tab 0.006 \tab 1.008 \cr
  0.0 \tab 2000 \tab 0.100 \tab 0.053 \tab 0.009 \tab 1.009 \cr
  0.0 \tab 4000 \tab 0.092 \tab 0.043 \tab 0.004 \tab 1.045 \cr
  0.1 \tab 500 \tab 0.422 \tab 0.298 \tab 0.116 \tab 1.027 \cr
  0.1 \tab 1000 \tab 0.629 \tab 0.492 \tab 0.270 \tab 0.991 \cr
  0.1 \tab 2000 \tab 0.860 \tab 0.787 \tab 0.578 \tab 0.983 \cr
  0.1 \tab 4000 \tab 0.984 \tab 0.969 \tab 0.897 \tab 0.970 \cr
  0.2 \tab 500 \tab 0.879 \tab 0.807 \tab 0.611 \tab 0.976 \cr
  0.2 \tab 1000 \tab 0.986 \tab 0.976 \tab 0.902 \tab 0.951 \cr
  0.2 \tab 2000 \tab 1.000 \tab 1.000 \tab 0.998 \tab 0.971 \cr
  0.2 \tab 4000 \tab 1.000 \tab 1.000 \tab 1.000 \tab 0.974
}

\bold{The test holds its size at every sample size tried}, and power rises monotonically to one at both alternatives. The last column is the diagnostic behind that: it is the mean \emph{reported} standard error of \eqn{\hat\rho} divided by the actual sampling standard deviation of \eqn{\hat\rho} across replications. The statistic is \eqn{\hat\rho^2/\widehat{\mathrm{Var}}(\hat\rho)}, so a rejection rate away from nominal comes either from \eqn{\hat\rho} itself or from the variance, and only the second would be a defect here. At 0.95 to 1.05 throughout, the variance estimate is sound -- which is what separates this test from \code{\link{uhet_test}}, where the same ratio is 0.93 and the test is conservative in consequence.

At the paper's own design (\eqn{N = 100}, \eqn{T = 20}) the larger alternatives it reports are saturated: power is 1.000 at every level at both \eqn{\rho = 0.4} and \eqn{\rho = 0.8}. The \eqn{\rho = 0.1} and \eqn{0.2} rows are included precisely because a column of ones says nothing about how power builds. Replication scripts: \code{spec_tests/hou2025_table2.R} and \code{spec_tests/hou2025_endog_size.R}.

\bold{Reading a non-rejection.} Failing to reject is not evidence for exogeneity. \eqn{\rho} is identified through the instruments, so weak instruments make \eqn{\hat\rho} imprecise and the test powerless; the print method says so rather than letting the p-value speak alone.
}
\references{
Hou, Z., Ramalho, J. J. S. and Roseta-Palma, C. (2025). Dealing with endogeneity in stochastic frontier models: A comparative assessment of estimators. \emph{Energy Economics} 151, 108922.

Amsler, C., Prokhorov, A. and Schmidt, P. (2016). Endogeneity in stochastic frontier models. \emph{Journal of Econometrics} 190(2), 280-288.

Kutlu, L. (2010). Battese-Coelli estimator with endogenous regressors. \emph{Economics Letters} 109(2), 79-81.
}
\seealso{\code{\link{ivsfm}} for the fit, \code{\link{inefficiency_test}} for the boundary-null case, \code{\link{spec_test}} and \code{\link{moment_range}} for the distributional assumptions.}
\examples{
set.seed(9)
n  <- 600
z  <- rnorm(n); x1 <- rnorm(n); eps <- rnorm(n)
x2 <- z + eps
v  <- 0.6 * eps + sqrt(1 - 0.36) * rnorm(n)
d  <- data.frame(y = 0.5 * x1 + 0.5 * x2 + v - abs(rnorm(n)),
                 x1 = x1, x2 = x2, z = z)

fit <- ivsfm(y ~ x1 + x2, endogenous = ~ x2, instruments = ~ z,
             data = d, model_name = "IVCF")
endogeneity_test(fit)
}
\keyword{htest}
