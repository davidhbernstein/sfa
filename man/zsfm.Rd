\name{zsfm}
\alias{zsfm}
\title{Zero-Inefficiency Stochastic Frontier Models}
\description{
Fits the zero-inefficiency stochastic frontier model of Kumbhakar, Parmeter and Tsionas (2013): the sample is a mixture of fully efficient firms (\eqn{u = 0}) and firms drawn from an ordinary normal-half-normal frontier.
\itemize{
\item \code{"ZISF"}: the probability of being fully efficient is one constant for the whole sample.
\item \code{"ZISF_Z"}: that probability depends on covariates, given in the second part of the formula (\code{y ~ x | z}).
}
}
\usage{zsfm(formula, model_name = c("ZISF", "ZISF_Z"), 
data, maxit.bobyqa = 10000,maxit.psoptim = 1000, maxit.optim = 1000, 
REPORT = 1, trace = 0, pgtol = 0,start_val = FALSE,PSopt = FALSE, 
optHessian = TRUE, inefdec = TRUE, upper = NA, 
Method = "L-BFGS-B",logit = TRUE,verbose=FALSE,
keep_objective = FALSE,rand.psoptim = NULL)}
\arguments{
\item{formula}{A symbolic description of the frontier, \code{y ~ x1 + x2}. For \code{"ZISF_Z"} a second part gives the covariates of the probability of full efficiency: \code{y ~ x1 + x2 | z1}.}
\item{model_name}{\code{"ZISF"} or \code{"ZISF_Z"}; see \sQuote{Description}. Matching ignores case, so \code{"zisf"} and \code{"ZISF"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine.}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine.}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine.}
\item{REPORT, trace, pgtol, upper}{Accepted but not used by \code{zsfm()}: the optimizer's own tracing is controlled by \code{verbose}, and its bounds are set internally.}
\item{start_val}{\code{FALSE} (the default) to build starting values from the data, or a numeric vector giving them directly, in the order of the rows of \code{out}.}
\item{PSopt}{Logical. Run the particle-swarm stage (\code{pso::psoptim}) between bobyqa and \code{optim}? Defaults to \code{FALSE}.}
\item{optHessian}{Logical. Finish with \code{optim()} and a numerically differentiated Hessian? Standard errors and t-values require \code{TRUE}, the default.}
\item{inefdec}{Logical. \code{TRUE} (the default) for a production frontier, \eqn{\varepsilon = v - u}; \code{FALSE} for a cost frontier, \eqn{\varepsilon = v + u}.}
\item{Method}{Method passed to \code{optim} for the final stage. Must be
\code{"L-BFGS-B"}: it is the only \code{optim} method that applies the lower
and upper bounds this package uses to keep scale parameters positive, and any
other value is refused with an error rather than quietly optimized without
them.}
\item{logit}{Logical, \code{"ZISF_Z"} only. The link from \eqn{z'\gamma} to the probability of full efficiency; see \sQuote{The probability of full efficiency}. Prefer the default, \code{TRUE}.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{keep_objective}{Logical. Retain the likelihood closure on the fitted object as \code{$objective}, so the score-based tools can difference it: \code{vcov(type = "bhhh")} (OPG/BHHH standard errors), \code{\link{influence_sfa}}, the \pkg{sandwich} methods, \code{\link{TIC}} and \code{vuong}. New in 1.2.1. \code{FALSE} by default, because a closure carries its enclosing environment, so a fit saved with one serializes the estimation data alongside the results.}
\item{rand.psoptim}{Integer. Seed for replication of psoptim. Defaults to \code{NULL}.}
}
\details{
\subsection{The model}{
With probability \eqn{p} a firm sits exactly on the frontier, so its composed error is noise alone, \eqn{\varepsilon = v}. With probability \eqn{1 - p} it is drawn from an ordinary normal-half-normal frontier. This is the natural specification when some fraction of the sample is thought to be fully efficient, which an ordinary frontier cannot represent: it forces \eqn{u > 0} almost surely.
}

\subsection{The probability of full efficiency}{
Under \code{"ZISF"}, \eqn{p} is one constant, parameterized as \eqn{p = \exp(-|\gamma|)}. That parameterization is exactly symmetric in \eqn{\gamma}: \eqn{+\gamma} and \eqn{-\gamma} give identical fits and the optimizer may return either, so the sign of the reported \code{gamma} carries no meaning.

Under \code{"ZISF_Z"}, \eqn{p} varies with the covariates \eqn{z} through the link chosen by \code{logit}:
\itemize{
\item \code{logit = TRUE} (the default): the logistic link \eqn{\Lambda(z'\gamma)}, which covers \eqn{(0,1)}.
\item \code{logit = FALSE}: \eqn{\Phi(z'\gamma)/(1 + \Phi(z'\gamma))}, which is \bold{bounded above by 0.5}. Under that link no observation can be assigned a probability of full efficiency above one half, whatever the data. It is not the probit analogue of the default and is retained only for backward compatibility.
}
}

\subsection{Efficiency predictions}{
Two per-observation predictions are returned, and they answer different questions:
\itemize{
\item \code{post.prob} is the posterior probability that the observation is fully efficient, given its residual.
\item \code{jlms} is the Jondrow et al. (1982) prediction \eqn{E[u \mid \varepsilon]} \emph{conditional on the observation being in the inefficient regime}. It is the ordinary normal-half-normal formula and does not weight by \code{post.prob}.
}
The prediction of \eqn{u} over both regimes is \code{(1 - post.prob) * jlms}. The two can differ by orders of magnitude when the fitted \eqn{p} is large, so use the product when you want \eqn{E[u \mid \varepsilon]} itself.
}

\subsection{Iteration defaults}{
Like \code{sfm()}, \code{zsfm()} has a closed-form likelihood, so its optimization defaults (thousands of iterations) are much higher than \code{psfm()}'s, which uses simulated maximum likelihood for several models. See \code{\link{psfm}}'s Details.
}
}
\references{Kumbhakar, S.C., Parmeter, C.F. and Tsionas, E.G. (2013) 'A zero inefficiency stochastic frontier model', \emph{Journal of Econometrics}, 172(1), pp. 66-76. doi:10.1016/j.jeconom.2012.08.021.

Jondrow, J., Lovell, C.A.K., Materov, I.S. and Schmidt, P. (1982) 'On the estimation of technical inefficiency in the stochastic frontier production function model', \emph{Journal of Econometrics}, 19(2-3), pp. 233-238.}
\author{Christopher F. Parmeter and David H. Bernstein}
\note{Standard errors require \code{optHessian = TRUE}.}
\seealso{\code{\link{sfm}}, \code{\link{psfm}}, \code{\link{ttsfm}}, \code{\link{sfareg-methods}}, \code{\link{panel89}}}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values. Rows are \code{gamma}, \code{sigv}, \code{sigu} and the frontier coefficients for \code{"ZISF"}; \code{sigv}, \code{sigu}, the frontier coefficients and then the \eqn{\gamma} coefficients for \code{"ZISF_Z"}.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{\code{"ZISF"} or \code{"ZISF_Z"}.}
\item{formula}{The formula used in the model specification.}
\item{jlms}{The Jondrow et al. (1982) inefficiency prediction \eqn{E[u \mid \varepsilon]} within the inefficient regime; see \sQuote{Efficiency predictions}.}
\item{post.prob}{Posterior probability that each observation is fully efficient.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{nobs}{The number of observations used in the fit.}
\item{objective}{The likelihood closure, when \code{keep_objective = TRUE}.}
}
\examples{
\donttest{
library(sfa)  

eqz     <- y ~ q1 + q2 + q3 + q4 + q5 + w1 + w2 + w3 + w4 | z

data(panel89)

f <- zsfm(formula    = eqz,
          model_name = "ZISF_Z",
          data       = panel89,
          logit      = TRUE)
f

## Inefficiency over both regimes, not just the inefficient one
head((1 - f$post.prob) * f$jlms)
}
}
