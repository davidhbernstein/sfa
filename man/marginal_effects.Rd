\name{marginal_effects}
\alias{marginal_effects}
\title{Marginal Effects of the Variance Determinants on Inefficiency}
\description{
For a model whose inefficiency scale depends on covariates (the \code{_Z}
models), the marginal effects of those covariates on the mean and variance of
the inefficiency term, \eqn{\partial E[u]/\partial z_k} and
\eqn{\partial \mathrm{Var}[u]/\partial z_k}.
}
\usage{
marginal_effects(object, average = FALSE)
}
\arguments{
  \item{object}{an object of class \code{"sfareg"} from \code{\link{sfm}},
    fitted with a variance-determinant formula such as
    \code{y ~ x1 + x2 | z}.}
  \item{average}{if \code{TRUE}, return only the named vector of average
    marginal effects. Default \code{FALSE}, which returns the full
    per-observation table with that vector attached as the \code{"average"}
    attribute.}
}
\details{
A \eqn{\delta} coefficient is not interpretable on its own: it sits in a log
link for a scale parameter, so its units are neither those of \eqn{u} nor
those of \eqn{z}, and its magnitude is not comparable across models whose link
differs. The marginal effect is on the scale of \eqn{u}, and is what applied
work reports.

Writing \eqn{s = \sigma_u(z)}, every derivative reduces to
\eqn{\partial s/\partial z_k} times a constant that cancels back into the
moment itself. For the half-normal and the exponential alike:

\deqn{\textrm{SD link:} \quad \partial E[u]/\partial z_k = \delta_k E[u],
      \qquad \partial \mathrm{Var}[u]/\partial z_k = 2\delta_k \mathrm{Var}[u]}
\deqn{\textrm{variance link:} \quad \partial E[u]/\partial z_k = (\delta_k/2) E[u],
      \qquad \partial \mathrm{Var}[u]/\partial z_k = \delta_k \mathrm{Var}[u]}

\strong{The factor of two matters.} \code{sfm()}'s \code{"NHN_Z"} and
\code{"NE_Z"} place the linear predictor on the standard deviation,
\eqn{\sigma_u = \exp(z'\delta)}, while \code{psfm()}'s \code{"TRE_Z"} and
\code{"GTRE_Z"} place it on the variance,
\eqn{\sigma_u = \sqrt{\exp(z'\delta)}}. Reading a \eqn{\delta} from one family
with the other's convention in mind is wrong by exactly this factor. Reporting
the marginal effect instead of the coefficient removes the trap, since the
effect is on the scale of \eqn{u} either way.

Effects are per observation, because \eqn{\sigma_u} varies with \eqn{z}.
Constant columns of the \eqn{z} design are dropped: the derivative with
respect to an intercept is not a marginal effect.

Standard errors are deliberately not reported. The effect is a nonlinear
function of both \eqn{\delta} and \eqn{z}, so its sampling distribution
requires the delta method through the full covariance matrix or a bootstrap;
a number printed without one would be read as inference.
}
\value{
With \code{average = FALSE}, a data frame with one row per observation:
\code{sigma_u}, \code{E_u}, \code{Var_u}, then one \code{dE_u.d<z>} and one
\code{dVar_u.d<z>} column per non-constant variance determinant. Attributes
\code{"average"} (the average marginal effects), \code{"link"} and
\code{"family"} record the convention used.

With \code{average = TRUE}, just the named vector of average marginal effects.
}
\seealso{
\code{\link{sfm}}, \code{\link{efficiency_ci}}
}
\examples{
dat <- data_gen_cs(N = 300, rand = 1, sig_u = 1, sig_v = 0.3,
                   cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
fit <- sfm(y_pcs_z ~ x1 + x2 | z, data = dat, model_name = "NHN_Z")

me <- marginal_effects(fit)
head(me)

## The conventional summary.
marginal_effects(fit, average = TRUE)
}
\keyword{models}
