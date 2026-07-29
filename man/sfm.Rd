\name{sfm}
\alias{sfm}
\title{sfm}
\description{Implementation of the cross-sectional stochastic frontier model across an array of distributional assumptions for both v and u (user specified). For panel models, see the psfm() call.}
\usage{sfm(formula, model_name, data,maxit.bobyqa,maxit.psoptim,maxit.optim,REPORT,
trace,pgtol,start_val,PSopt,optHessian,inefdec,upper,Method,
robust=c("mle","mlqe","psi","mdpd"),c_mlqe=0.20,eta,alpha,verbose=FALSE,
rand.psoptim=NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{model name for the estimation includes the: normal-half normal (NHN),
normal-exponential (NE), student's t-half t (THT),  Normal-Rayleigh (NR),
and the normal-truncated normal (NTN).}
\item{data}{A data set}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{trace}
\item{pgtol}{pgtol}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T or F)}
\item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
\item{inefdec}{Production or cost function}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See 'Details' within optim.}
\item{robust}{One of \code{"mle"} (default), \code{"mlqe"}, \code{"psi"}, or
\code{"mdpd"}. Selects a robust divergence-based objective in place of
ordinary maximum likelihood -- see Details. Currently only implemented for
\code{model_name = "NHN"}; \code{sfm()} stops with a clear error if a
non-\code{"mle"} value is requested for any other model.}
\item{c_mlqe}{Numeric, default \code{0.20}. The MLqE power parameter
(\code{c} in \code{f^c}; the corresponding \code{q = 1-c}), used only when
\code{robust = "mlqe"}.}
\item{eta}{Numeric, default \code{0.01}. The Psi-likelihood power parameter,
used only when \code{robust = "psi"}.}
\item{alpha}{Numeric, default \code{0.2}. The MDPD (minimum density power
divergence) power parameter, used only when \code{robust = "mdpd"}.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
}
\details{The options include the Normal-Half Normal (NHN), Normal-exponential (NE), Student's t-Half t (THT),  and the Normal-Truncated Normal (NTN).  NHN_Z and NE_Z are extensions for the NHN and NE models that allow for modeling the u-component of those models with determinants of inefficiency.

Outputs include E[exp(-u)|e] given by exp_u_hat, following Battese and Coelli (1988, JoE), where appropriate.

Note on default iteration counts: \code{sfm()}'s optimization defaults (thousands of iterations) are much higher than \code{psfm()}'s (100/10/10). This is intentional: \code{sfm()}'s likelihoods are closed-form (no simulation), so each objective-function evaluation is cheap and many default iterations cost little; several of \code{psfm()}'s models require simulated maximum likelihood over Halton draws, where each evaluation is far more expensive. See \code{\link{psfm}}'s Details for more.

\strong{Robust divergence-based estimation (\code{robust} argument).} For
\code{model_name = "NHN"}, \code{sfm()} can maximize a robust divergence
objective instead of the ordinary log-likelihood, downweighting the
influence of outlying observations:
\itemize{
  \item \code{"mlqe"}: minimum L\eqn{q} estimation (Lq, Ferrari and Yang
  2010), power parameter \code{c_mlqe} (\eqn{q = 1-c}). Objective
  \eqn{\sum (f(e_i)^c - 1)/c}.
  \item \code{"psi"}: Psi-likelihood (Basu et al.-type divergence family),
  power parameter \code{eta}. Objective
  \eqn{\sum f(e_i)^c/c - n \, I_{1+c}/(1+c)}, where
  \eqn{I_{1+c} = \int f(e;\theta)^{1+c}\,de} is evaluated numerically at the
  current parameter values on every optimizer iteration.
  \item \code{"mdpd"}: minimum density power divergence (Basu et al. 1998),
  power parameter \code{alpha}. Objective
  \eqn{((1+c)/c) \sum f(e_i)^c - n \, I_{1+c}}.
}
All three reduce to \code{"mle"} as their power parameter approaches 0.
\strong{Psi and MDPD produce identical parameter point estimates for the
same power parameter} (\eqn{c = \eta = \alpha}) -- the two objectives are a
positive rescaling of one another, \eqn{\text{MDPD}(c) = (1+c)\cdot
\text{Psi}(c)}, which does not change the argmax; both are exposed for
compatibility with the corresponding paper's three-method terminology, not
because they differ numerically. MLqE is a genuinely distinct estimator.

Because these are M-estimators rather than ordinary MLE, the naive
Hessian-inverse standard error (\code{optHessian = TRUE}'s usual output) is
not statistically valid for them -- it assumes the information-matrix
equality, which these objectives do not satisfy in general. When
\code{robust != "mle"}, \code{sfm()} instead reports a sandwich-form
standard error, \eqn{A^{-1} B A^{-1}}, where \eqn{A} is the Hessian of the
summed robust objective (as already computed for the point estimate) and
\eqn{B = \sum_i g_i g_i'} is built from the per-observation gradient of that
objective (obtained via numerical differentiation,
\code{\link[numDeriv]{jacobian}}). This requires \code{optHessian = TRUE};
\code{std.errors}/\code{t.values} are \code{NA} otherwise, same as for
ordinary MLE.

Formulas verified against an independently audited reference implementation
(cross-checked by direct numerical integration and finite differences),
including the sandwich SE's own verification (the c -> 0 limit against the
ordinary MLE Hessian SE, and a Monte Carlo check against the empirical
spread of repeated fits). Currently only implemented for
\code{model_name = "NHN"} -- other models will error clearly if
\code{robust != "mle"} is requested.}
\author{David H. Bernstein and Alexander Stead}
\note{Standard errors require optHessian set to TRUE. When
\code{robust != "mle"} (\code{model_name = "NHN"} only), a sandwich-form SE
is reported instead of the ordinary Hessian-inverse one -- see Details.}

\seealso{see also}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{exp_u_hat}{Predicted technical efficiency (expected values). Available for models: NHN, NHN_Z, NR, NG, and NNAK.}
\item{med_u_hat}{Predicted technical efficiency (median values). Available only for the NHN model.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}). Sandwich-form (not the ordinary Hessian-inverse form) when \code{robust != "mle"} -- see Details.}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{robust}{Only present for \code{model_name = "NHN"}: the \code{robust} argument used to fit this model (\code{"mle"}, \code{"mlqe"}, \code{"psi"}, or \code{"mdpd"}).}
\item{robust_c}{Only present for \code{model_name = "NHN"}: the power parameter (\code{c_mlqe}/\code{eta}/\code{alpha}, whichever applies) used, or \code{NA} when \code{robust = "mle"}.}
}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN_Z",
                    data       = cs_data_trial,          PSopt      = TRUE)

## Robust divergence-based estimation (NHN only): downweights the influence
## of outlying observations relative to ordinary MLE. std.errors/t.values
## are NA for these methods -- see Details.
cs.nhn.mlqe <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NHN",
                    data = cs_data_trial, robust = "mlqe", c_mlqe = 0.20)
cs.nhn.psi  <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NHN",
                    data = cs_data_trial, robust = "psi", eta = 0.20)
}
}