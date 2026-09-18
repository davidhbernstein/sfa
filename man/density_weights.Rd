\name{density_weights}
\alias{density_weights}
\title{Density-power weights from a robust frontier fit}
\description{
Returns the weight each observation receives in the estimating equation of a
maximum \eqn{L_q}-likelihood or minimum density-power divergence fit. The weight
is the fitted density raised to the power \eqn{c}, normalised so the largest is
one. It is formed as \eqn{\exp(c \log f)} rather than \eqn{(\exp \log f)^c},
so an observation whose density underflows the natural scale still receives its
own weight and keeps its place in the ordering.
}
\usage{
density_weights(object, sigma_v = NULL, sigma_u = NULL, c = NULL,
                normalize = TRUE, log = FALSE)
}
\arguments{
\item{object}{A fitted frontier model, or a numeric vector of composed
residuals.}
\item{sigma_v, sigma_u}{Scale parameters. Required when \code{object} is a
numeric vector; taken from the fit otherwise.}
\item{c}{Robustness tuning parameter. Taken from the fit when available.}
\item{normalize}{Logical; scale so the maximum weight is one. Default
\code{TRUE}.}
\item{log}{Logical; return \eqn{\log w} instead of \eqn{w}. Default
\code{FALSE}. Use this to rank or threshold observations: the weights of the
most heavily discounted observations span hundreds of orders of magnitude and
are not all distinguishable once exponentiated.}
}
\value{A numeric vector of weights in \eqn{(0, 1]}, one per observation, or of
log weights in \eqn{(-\infty, 0]} when \code{log = TRUE}.}
\section{What this diagnostic can and cannot find}{
The weight is a function of the fitted density at the observation, so it
responds to observations lying far from the fitted surface. It has no purchase
on an observation that is wrong in a \emph{regressor} and, because the fitted
surface bends toward it, ends up close to that surface. A mis-recorded input can
therefore be highly influential and still receive a weight near one. This is not
a defect of the implementation: the bounded influence of these estimators is
bounded with respect to \emph{response} contamination, conditional on a fixed or
bounded design. Pair this with a leverage diagnostic, and with case deletion,
before concluding that an observation is or is not driving a result.
}
\examples{
set.seed(1)
e <- rnorm(200, 0, 0.3) - abs(rnorm(200, 0, 0.6))
w <- density_weights(e, sigma_v = 0.3, sigma_u = 0.6, c = 0.217)
summary(w)
head(order(w), 10)   # the observations the estimator discounts most

## Rank on the log scale when the fit is sharply peaked: there the smallest
## weights are far below the smallest positive double and cannot be told
## apart once exponentiated.
lw <- density_weights(e, sigma_v = 0.02, sigma_u = 0.9, c = 0.217, log = TRUE)
head(order(lw), 10)
}
\references{
Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of the Stochastic Frontier Model. Working paper.
}
\seealso{\code{\link{hscore_select}}, \code{\link{calibrate_c}}}
