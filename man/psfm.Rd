\name{psfm}
\alias{psfm}
\title{psfm}
\description{Function to implement various panel data stochastic frontier estimators}
\usage{psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE",
                    "GTRE", "TFE", "FD", "GTRE_SEQ1", "GTRE_SEQ2", "SSFE",
                    "PL80", "BC92"), data,
                    maxit.bobyqa = 100, maxit.psoptim = 10, maxit.optim =
                    10, REPORT = 1, trace = 3, pgtol = 0, individual,
                    halton_num = NULL, start_val = FALSE, gamma = FALSE,
                    PSopt = FALSE, optHessian = TRUE, inefdec = TRUE, Method = "L-BFGS-B",
                    verbose = FALSE, rand.gtre = NULL, rand.psoptim = NULL,
                    OPG_calc = FALSE, time = NULL)}
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{model name for the estimation}
  \item{data}{a pdata.frame}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
  \item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
  \item{REPORT}{reporting parameter}
  \item{trace}{trace}
  \item{pgtol}{pgtol}
  \item{individual}{individual unit in the regression model}
  \item{halton_num}{number of Halton draws to use in SML models}
  \item{start_val}{starting value (optional)}
  \item{gamma}{gamma}
  \item{PSopt}{use psoptim optimization routine (T or F)}
  \item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
  \item{inefdec}{Production or cost function}
  \item{Method}{The method to be used for optim. See 'Details' within optim.}
  \item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
  \item{rand.psoptim}{Integer. Seed for replication of psoptim.  Default to \code{NULL}.}
  \item{rand.gtre}{Integer. Seed for replication of the gtre model.  Default to \code{NULL}.}
  \item{OPG_calc}{Logical. Should the outer product of gradients be calculated for inference.  Default to \code{FALSE}.}
  \item{time}{Character. Name of a time-period column in \code{data}, used only by the \code{"PL80"}/\code{"BC92"} models. If \code{NULL} (the default), a sequential within-individual time index is constructed from row order, which assumes \code{data} is already sorted chronologically within each individual.}
}
\details{The generalized true random effects model (GTRE, 4-component model) and true random effects models (TRE)  are both estimated by simulated maximum likelihood based on the paper by the Fillipini and Greene (2016, JPA). The TRE_Z and GTRE_Z allow for modeling the u-component of the GTRE and TRE with determinants of inefficiency. The first-difference estimator (FD) of Wang and Ho (2010, JoE) as well as the True Fixed Effect model estimated by within-maximum likelihood of Chen, Schmidt and Wang (2014, JoE) are also available.

\code{"SSFE"} is the classical Schmidt and Sickles (1984, JBES) fixed-effects (LSDV) estimator: a standard within/fixed-effects panel regression, with technical inefficiency read directly off the estimated firm effects (each firm's distance from the best-performing firm in the sample) rather than estimated by maximum likelihood. It requires no numerical optimization and returns no log-likelihood -- \code{logLik()}/\code{AIC()}/\code{BIC()} are not defined for it (see \code{\link[sfa]{sfareg-methods}}). Note this is a different model from \code{"TFE"} above, despite the similar name: \code{"TFE"} is Chen, Schmidt and Wang's (2014) maximum-likelihood true fixed-effects estimator, which allows a separate stochastic noise term on top of the fixed effects, while \code{"SSFE"} attributes all within-firm deviation to the two-sided error/inefficiency directly.

\code{"PL80"} (time-invariant panel inefficiency, Pitt and Lee, 1980, JoE) and \code{"BC92"} (time-varying decay, Battese and Coelli, 1992, JPA) are both natively estimated maximum-likelihood implementations of the same "error components frontier": y_it = x_it'beta + v_it - B_it*u_i (production; sign flips for cost), with v_it iid normal, u_i iid half-normal and constant across time within a firm, and B_it = 1 for \code{"PL80"} or B_it = exp(-eta*(t - Tref)) for \code{"BC92"} (Tref being the last time period in the whole panel). The closed-form log-likelihood is obtained by integrating u_i out of the joint density of (epsilon_i, u_i); it nests the ordinary cross-sectional normal-half-normal likelihood when every firm has a single observation. Verified against \code{frontier::sfa()} (matching coefficients, log-likelihood, and \code{frontier::efficiencies()}-style predicted efficiencies to numerical precision across balanced/unbalanced panels and production/cost specifications) before that dependency was removed. The \code{time} argument above supplies each observation's period for \code{"BC92"}'s decay term; if left \code{NULL}, a sequential within-individual time index is constructed from row order, which assumes \code{data} is already sorted chronologically within each individual.

Note on default iteration counts: \code{psfm()}'s default \code{maxit.bobyqa}/\code{maxit.psoptim}/\code{maxit.optim} (100/10/10) are much lower than the defaults in \code{sfm()}/\code{zsfm()}/\code{ttsfm()} (thousands to tens of thousands). This is intentional, not an oversight: several \code{psfm()} models (GTRE, GTRE_Z, TRE, TRE_Z) are estimated by simulated maximum likelihood over Halton draws, making each objective-function evaluation far more expensive than the closed-form likelihoods used elsewhere in this package. If a \code{psfm()} fit does not appear to have converged, raise these values explicitly -- the low defaults trade some convergence thoroughness for a fast default call.}
\references{Fillipini and Greene (2016, JPA); Wang and Ho (2010, JoE); Chen, Schmidt and Wang (2014, JoE); Pitt and Lee (1980, JoE); Battese and Coelli (1992, JPA); Battese and Coelli (1988, JoE)}
\author{David Bernstein}
\note{Standard errors require optHessian set to TRUE}

\seealso{see also}
\value{
An object of class \code{"sfareg"} containing components that vary by model. All models return:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{model_name}{The name of the panel stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{data}{The data used in estimation.}

\bold{Additional model-specific components:}

\bold{For GTRE and GTRE_Z models:}
\item{H}{Predicted time-invariant technical efficiency for each individual.}

\bold{For GTRE_Z model:}
\item{out_opg}{OPG (outer product of gradients) standard errors and t-values along with parameters.}
\item{out_sandwich}{sandwich standard errors and t-values along with parameters.}

\bold{For GTRE, GTRE_Z, TRE and TRE_Z models:}
\item{U}{Predicted time-varying technical efficiency for each observation.}

\bold{For TFE model:}
\item{r_hat_m}{Estimated individual-specific random effects.}
\item{exp_u_hat}{Predicted technical efficiency.}

\bold{For FD model:}
\item{u_hat}{Predicted technical efficiency in levels.}
\item{h_hat}{Estimated z heterogeneity function values.}
\item{exp_u_hat}{Predicted technical efficiency.}

\bold{For GTRE_SEQ1 and GTRE_SEQ2 models:}
\item{other_parms}{A matrix of additional parameters (lambda, sigma, beta_0 for SEQ1; sigma_u, sigma_v, sigma_h, sigma_r, lambda, sigma for SEQ2).}

\bold{For SSFE model:}
\item{alpha_hat}{The estimated firm-specific fixed effects (named by individual).}
\item{u_hat}{Estimated technical inefficiency for each firm, as its distance from the best-performing firm's fixed effect.}
\item{exp_u_hat}{Predicted technical efficiency, \code{exp(-u_hat)}.}

\bold{For PL80 and BC92 models:}
\item{exp_u_hat}{Predicted technical efficiency for each observation, E[exp(-B_it*u_i)|epsilon_i] (generalizes Battese and Coelli's (1988, JoE) point estimator to the panel case, pooling all of a firm's observations to predict u_i and scaling by B_it per period).}
}

\note{
The GTRE_SEQ1 and GTRE_SEQ2 models use sequential estimation methods and do not return optimization objects or starting values. The SSFE model is likewise not maximum-likelihood based and returns no optimization object, starting values, or log-likelihood. All panel models require the \code{individual} argument to identify panel units.
}
\examples{
\donttest{
library(sfa)     

data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                         sig_u = 1,  sig_v = 0.3, 
                         sig_r = .2, sig_h = .4, 
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

max_tre_z   <-  psfm(formula    = y_tre_z ~ x1 +x2| z_gtre,
                     model_name = "TRE_Z",
                     data       = data_trial,
                     individual = "name",
                     PSopt      = TRUE)
}
}