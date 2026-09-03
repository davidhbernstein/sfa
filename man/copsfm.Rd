\name{copsfm}
\alias{copsfm}
\title{Stochastic Frontier with Dependence Between Noise and Inefficiency}
\description{
Fits a cross-sectional stochastic frontier in which the noise \eqn{v} and the
inefficiency \eqn{u} are \emph{dependent}, with the dependence carried by a
copula. Every other model in the package assumes independence.
}
\usage{
copsfm(formula, data, copula = c("gaussian", "fgm"), inefdec = TRUE,
       n_nodes = 128, maxit.bobyqa = 10000, maxit.psoptim = 1000,
       maxit.optim = 1000, start_val = FALSE, PSopt = FALSE,
       optHessian = TRUE, Method = "L-BFGS-B", verbose = FALSE,
       rand.psoptim = NULL)
}
\arguments{
\item{formula}{Two-sided formula for the frontier, \code{y ~ x1 + x2}. A
  \code{|} segment is an error: variance determinants are not yet supported
  alongside a copula.}
\item{data}{A \code{data.frame}.}
\item{copula}{\code{"gaussian"} (parameter \code{rho}) or \code{"fgm"}
  (Farlie--Gumbel--Morgenstern, parameter \code{theta}). Both reduce to
  independence at 0.}
\item{inefdec}{\code{TRUE} for a production frontier, \code{FALSE} for cost.}
\item{n_nodes}{Gauss--Legendre nodes for the integral over \eqn{u}.}
\item{maxit.bobyqa, maxit.psoptim, maxit.optim, start_val, PSopt, optHessian,
  Method, verbose, rand.psoptim}{As in \code{\link{sfm}}.}
}
\details{
The package's taxonomy already writes the general joint density as
\deqn{f_{V,U}(v,u) = f_V(v)\,f_U(u)\,c(F_V(v), F_U(u); \rho),}
with every existing specification setting \eqn{c \equiv 1}. This function
relaxes exactly that and nothing else: the marginals remain normal and
half-normal. Composing over \eqn{v = \varepsilon + Su},
\deqn{f_\varepsilon(\varepsilon) = \int_0^\infty f_V(\varepsilon + Su)\,f_U(u)\,
      c\!\left(F_V(\varepsilon + Su), F_U(u)\right) du,}
a one-dimensional integral evaluated by Gauss--Legendre quadrature on
\eqn{u = t/(1-t)} rather than by simulation.

\bold{Why \code{n_nodes} defaults to 128.} Measured against the closed-form
normal/half-normal density at independence, the largest absolute error in
\eqn{\log f} is 1.3e-2 at 32 nodes, 2.5e-5 at 64, and 8.2e-14 at 128. Sixty-four
looks adequate and is not: 1e-5 per observation is 0.02 in a log-likelihood over
2000 points, which is the scale at which competing modes are compared.

\bold{The dependence parameter needs a large sample. This is the thing to know
before using it.} The frontier slopes behave normally, but \eqn{\rho} does not.
On a Gaussian-copula design with true \eqn{\rho = 0.6}, \eqn{\sigma_U = 1},
\eqn{\sigma_V = 0.4}:

\tabular{lrrr}{
\eqn{n} \tab \eqn{\hat\rho} \tab \eqn{\hat\sigma_U} \tab \eqn{\hat\sigma_V} \cr
2000    \tab -0.247 \tab 0.571 \tab 0.275 \cr
8000    \tab  0.528 \tab 0.921 \tab 0.379 \cr
20000   \tab  0.571 \tab 0.971 \tab 0.389
}

It is consistent -- the estimates converge on the truth -- but at \eqn{n = 2000}
the sign can come out wrong. The cause is joint identification rather than a
bad optimiser: with \eqn{\beta} and the scales held at their true values the
likelihood peaks exactly at \eqn{\rho = 0.6}, so \eqn{\rho} is well identified
\emph{conditionally}; it is the trade-off against \eqn{\sigma_U}, \eqn{\sigma_V}
and the intercept that defeats it in moderate samples, and the fitted point
genuinely attains a higher likelihood than the truth.

Treat \eqn{\hat\rho} as informative only in large samples, and compare against
an independent fit (\code{\link{sfm}}) before reading anything into it.

\bold{And below a few hundred observations, the independent fit is often the better
estimate -- even when the dependence is real.} Fitting the true copula gives a nearly
unbiased but very noisy \eqn{\hat\sigma_U}; ignoring the copula gives a precise but
biased one. RMSE of \eqn{\hat\sigma_U} (truth 1.0) on Gaussian-copula data, 25
replications:

\tabular{lrrr}{
design \tab \eqn{n} \tab true copula \tab independence \cr
\eqn{\rho = 0.5}  \tab 400  \tab 0.431 \tab \bold{0.240} \cr
\eqn{\rho = 0.5}  \tab 1000 \tab \bold{0.151} \tab 0.220 \cr
\eqn{\rho = -0.5} \tab 400  \tab 0.251 \tab \bold{0.196} \cr
\eqn{\rho = -0.5} \tab 1000 \tab \bold{0.175} \tab 0.203
}

The crossover sits between \eqn{n = 400} and \eqn{n = 1000}. Below it the variance
cost of estimating \eqn{\rho} jointly with \eqn{\sigma_U} exceeds the bias it
removes, so plain \code{\link{sfm}} wins on RMSE; above it \code{copsfm} wins. The
independence bias does not shrink with \eqn{n} (about 22\% of \eqn{\sigma_U}, since
it is misspecification rather than noise), while the copula fit's RMSE falls as it
should.

\bold{Node counts, measured in paired runs.} Forty datasets sent to every node count,
so the differences are the quadrature and not the draw. Against the 128-node default,
the mean difference in \eqn{\hat\rho} is +0.0002 at 256 nodes (\eqn{t = 1.37}) and
+0.0046 at 64 (\eqn{t = 1.47}) -- neither distinguishable -- but \bold{-0.0783 at 32}
(\eqn{t = -2.79}, \eqn{p = 0.008}), with per-dataset differences as large as 0.52. So
the quadrature has converged by 64, and 32 must not be used as a fast option: it is
biased, not merely noisier. For scale, \eqn{\hat\rho}'s own sampling standard
deviation is 0.385 at \eqn{n = 600}, which dwarfs every node-count effect here.

\bold{Families.} Only densities that could be verified are offered. Each is
checked in the tests two ways: it integrates to 1 over the unit square, and it
returns exactly 1 at the independence parameter. Frank, Clayton and Gumbel are
deliberately absent rather than transcribed without a source to check against.
}
\value{
An object of class \code{"sfareg"}, with \code{copula} and \code{copula_par}
alongside the usual components. \code{jlms} is \eqn{E[u \mid \varepsilon]}
computed by the same quadrature the likelihood used.
}
\references{
Smith, M. D. (2008). Stochastic frontier models with dependent error
components. \emph{The Econometrics Journal}, \bold{11}(1), 172--192.
}
\seealso{\code{\link{sfm}} for the independent case.}
\examples{
\donttest{
set.seed(4)
n   <- 4000
x1  <- rnorm(n); x2 <- rnorm(n)
z1  <- rnorm(n); z2 <- 0.6 * z1 + sqrt(1 - 0.6^2) * rnorm(n)
y   <- 0.5 + 0.8 * x1 - 0.4 * x2 + 0.4 * z1 - qnorm((1 + pnorm(z2)) / 2)
dat <- data.frame(y = y, x1 = x1, x2 = x2)

fit <- copsfm(y ~ x1 + x2, data = dat)
fit$out
}
}
\keyword{models}
