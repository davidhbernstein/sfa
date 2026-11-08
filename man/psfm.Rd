\name{psfm}
\alias{psfm}
\title{Panel Stochastic Frontier Models}
\description{Function to implement various panel data stochastic frontier estimators}
\usage{psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE",
                    "GTRE", "GTRE_FML", "TFE", "TFE_WMLE", "FD",
                    "GTRE_SEQ1", "GTRE_SEQ2", "SSFE",
                    "PL80", "PL80_MVTN", "BC92", "K1990", "K1990modified",
                    "CSS", "LS", "SSRE", "SSCRE", "KSS"), data,
                    maxit.bobyqa = 5000, maxit.nlminb = 500,
                    maxit.psoptim = 100, maxit.optim = 1000,
                    REPORT = 1, trace = 3, pgtol = 0, individual,
                    halton_num = NULL, start_val = FALSE, gamma = FALSE,
                    PSopt = FALSE, optHessian = TRUE, inefdec = TRUE, Method = "L-BFGS-B",
                    verbose = FALSE, rand.gtre = NULL, rand.psoptim = NULL,
                    OPG_calc = FALSE,
                    estimator = c("fiml", "sml", "seq1", "seq2"),
                    collinear_action = c("start_only", "error", "warn_drop"),
                    time = NULL, tfe_lambda_max = 100,
                    keep_objective = FALSE, kss_L = "auto",
                    kss_smooth = "auto", kss_L_max = 7L, mundlak = NULL,
                    sim_type = c("halton", "sobol", "torus", "uniform"),
                    antithetics = FALSE, sim_burn = NULL,
                    sim_scrambling = 0L, sim_prime = NULL)}
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{Which panel estimator to fit. All are described in \sQuote{Details}; in brief:
\describe{
  \item{\code{"TRE"}, \code{"TRE_Z"}}{Greene's (2005) true random effects, by simulated maximum likelihood. \code{_Z} adds determinants of \eqn{\sigma_u}.}
  \item{\code{"GTRE"}, \code{"GTRE_Z"}}{the four-component model separating persistent from transient inefficiency, by simulated maximum likelihood. \code{_Z} adds determinants.}
  \item{\code{"GTRE_FML"}}{the same four-component model by full information maximum likelihood through its closed-skew-normal representation; deterministic, but requires a balanced panel.}
  \item{\code{"GTRE_SEQ1"}, \code{"GTRE_SEQ2"}}{two-step/moment-based estimators of the four-component model. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"TFE"}}{Greene's (2005) true fixed effects. \bold{Changed meaning in 1.1.3}; see \sQuote{Details}.}
  \item{\code{"TFE_WMLE"}}{Chen, Schmidt and Wang's (2014) within maximum-likelihood estimator.}
  \item{\code{"FD"}}{the first-difference estimator of Wang and Ho (2010).}
  \item{\code{"SSFE"}}{Schmidt and Sickles' (1984) fixed-effects estimator, inefficiency read off the firm effects. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"CSS"}}{Cornwell, Schmidt and Sickles (1990): a firm-specific quadratic in time, so inefficiency varies over time and firms may cross over. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"LS"}}{Lee and Schmidt (1993): one common temporal pattern scaled per firm, \eqn{\alpha_{it} = \delta_t \alpha_i}. Time-varying but no crossover. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"SSRE"}}{Schmidt and Sickles (1984) random effects: the GLS estimator, with inefficiency read off the firm effects' BLUPs. More efficient than \code{"SSFE"} and identifies time-invariant regressors, but assumes the effects are uncorrelated with the regressors. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"SSCRE"}}{Correlated random effects: \code{"SSRE"} augmented with the within-firm means of the time-varying regressors (Mundlak 1978), so the correlation between effects and regressors is modelled. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"KSS"}}{Kneip, Sickles and Song (2012): the firm effect is a smooth function of time in an \eqn{L}-dimensional space whose basis is estimated from the data. Nests \code{"SSFE"}, \code{"LS"} and \code{"CSS"}. Balanced panels only. Not maximum likelihood: no \code{logLik()}.}
  \item{\code{"PL80"}}{Pitt and Lee (1980): time-invariant inefficiency.}
  \item{\code{"BC92"}}{Battese and Coelli (1992): exponential time decay.}
  \item{\code{"PL80_MVTN"}}{Pitt and Lee's (1981) Model III: inefficiency varies over time and is CORRELATED within a firm, through a multivariate normal truncated to the negative orthant. This is the likelihood from their Appendix 2, which they derived but did not use, calling the T-dimensional integrals intractable. Balanced panels only. See \code{\link{PL80_MVTN}}.}
  \item{\code{"K1990"}, \code{"K1990modified"}}{Kumbhakar (1990) time patterns for inefficiency.}
} Matching ignores case, so \code{"tre"}, \code{"TRE"} and \code{"Tre"} are the same choice.}
  \item{data}{A panel data set. This may be an ordinary \code{data.frame} (or tibble/data.table), in which case the panel index is built from the \code{individual} and \code{time} arguments, or a \code{plm::pdata.frame} that already carries an index, which is used as-is.}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.nlminb}{Maximum number of iterations for the nlminb optimization routine.}
  \item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
  \item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
  \item{REPORT, trace, pgtol}{Accepted but not used by \code{psfm()}: the optimizer's own tracing is controlled by \code{verbose}.}
  \item{individual}{individual unit in the regression model}
  \item{halton_num}{Number of simulation draws per firm for the simulated-ML models (\code{"GTRE"} with \code{estimator = "sml"}, \code{"TRE"}, \code{"GTRE_Z"}, \code{"TRE_Z"}). Despite the name, it counts draws of whichever sequence \code{sim_type} selects. \code{NULL} uses \eqn{\lceil\sqrt{n}\rceil + 100}.}
  \item{start_val}{starting value (optional)}
  \item{gamma}{Logical. Report the variance parameters in the \eqn{\gamma = \sigma_u^2/(\sigma_u^2 + \sigma_v^2)} parameterization instead of the default. Defaults to \code{FALSE}.}
  \item{PSopt}{use psoptim optimization routine (T or F)}
  \item{optHessian}{Logical. Finish with \code{optim()} and a numerically differentiated Hessian? Standard errors and t-values require \code{TRUE}, the default. With \code{FALSE} the final \code{optim} stage is not run at all, so the estimates reported are the preceding stage's and can differ from the default fit's wherever that stage had not already converged. \code{"PL80"}, \code{"BC92"}, \code{"K1990"} and \code{"K1990modified"} are the exception: their final stage runs either way, so \code{FALSE} costs only the Hessian.}
  \item{inefdec}{Production or cost function}
  \item{Method}{Method passed to \code{optim} for the final stage. Must be
  \code{"L-BFGS-B"}: it is the only \code{optim} method that applies the lower
  and upper bounds this package uses to keep scale parameters positive. Any
  other value is an error.}
  \item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
  \item{rand.psoptim}{Integer. Seed for the particle-swarm stage, so that a fit with \code{PSopt = TRUE} can be reproduced. It has no effect unless \code{PSopt = TRUE}. Varying it with \code{PSopt = TRUE} checks sensitivity to the starting point; see \code{rand.gtre}. Default \code{NULL}.}
  \item{rand.gtre}{Integer. Randomizes the simulation draws of the simulated-ML models listed under \code{halton_num}, by a uniform shift of the sequence (see \code{sim_type}). Fixing it makes a fit reproducible. Default \code{NULL}.

  \bold{Changing it changes the likelihood being maximized, not the starting point.} With the default \code{PSopt = FALSE} the start is deterministic, so refits under different \code{rand.gtre} start from the same place and maximize different simulated likelihoods. Their log-likelihoods are \emph{not comparable}. The refit with the highest one is not a better solution, only the one with the most favourable simulation error; at the default draw count such refits can even exceed the exact maximum that \code{estimator = "fiml"} gives.

  To measure this variation, use \code{\link{simulation_se}}, which refits under independent values of \code{rand.gtre} and reports the spread of the parameters and of the log-likelihood. To reduce it, raise \code{halton_num}. To check sensitivity to the \emph{starting point}, set \code{PSopt = TRUE} and vary \code{rand.psoptim}, or supply \code{start_val}. On a balanced panel, \code{model_name = "GTRE"} with \code{estimator = "fiml"} (the default) involves no simulation.}
  \item{OPG_calc}{Logical. Should the outer product of gradients be calculated for inference.  Default to \code{FALSE}.}
  \item{estimator}{How to estimate the four-component GTRE model, in the same spirit as \code{\link{sfm}}'s \code{estimator}: one model, several routes to it. Applies to \code{model_name = "GTRE"} only, and is ignored with a warning elsewhere.
\describe{
  \item{\code{"fiml"}}{\bold{the default.} Full information maximum likelihood through the closed-skew-normal representation. Deterministic: repeated fits to the same data agree exactly. Requires a \bold{balanced} panel.}
  \item{\code{"sml"}}{Simulated maximum likelihood over Halton draws. Handles unbalanced panels. This is what \code{model_name = "GTRE"} meant through version 1.1.3.}
  \item{\code{"seq1"}, \code{"seq2"}}{The two-step moment-based estimators. Not maximum likelihood, so they carry no \code{logLik()}.}
}
The older names remain and are unchanged: \code{"GTRE_FML"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} select the same three routes directly.

On an \bold{unbalanced} panel \code{"fiml"} cannot be fitted. If you took the default, \code{psfm()} warns and falls back to \code{"sml"}. If you asked for \code{"fiml"} explicitly it errors.}
  \item{collinear_action}{One of \code{"start_only"} (default), \code{"error"}, or \code{"warn_drop"}. Controls what happens when the \emph{between-individual} design matrix used for starting values is rank deficient. This is checked separately from the pooled design because the two can differ: time dummies, for example, are estimable pooled but collapse onto the intercept once averaged within each unit, which breaks the random-effects initialization. \code{"start_only"} keeps the model you asked for, drops the offending columns from the starting-value regression only, and takes their starts from pooled OLS; \code{"error"} stops and names the columns; \code{"warn_drop"} removes those terms from the model. All three report what was found. Only relevant to the models initialized by a random-effects regression (\code{GTRE}, \code{TRE} and their \code{_Z}/\code{SEQ} variants).}
  \item{time}{Character. Name of a time-period column in \code{data}, used by the error-components frontier models (\code{"PL80"}, \code{"BC92"}, \code{"K1990"}, \code{"K1990modified"}). If \code{NULL} (the default), a sequential within-individual time index is built from row order, which assumes \code{data} is sorted chronologically within each individual.}
  \item{tfe_lambda_max}{Numeric. Upper bound on \eqn{\lambda = \sigma_u/\sigma_v} when fitting \code{model_name = "TFE"}; ignored by every other model. Defaults to \code{100}, so \eqn{\sigma_v} cannot fall below about 1\% of \eqn{\sigma}. It guards against a degeneracy of Greene's likelihood (see \sQuote{Details}); it is not a convergence tolerance. A fit that pins at the bound warns and should not be read as an estimate of \eqn{\lambda}.}
\item{keep_objective}{Retain the likelihood closure on the fitted object as \code{$objective}, for every \code{model_name} fit by maximum likelihood: the simulated-ML models \code{"GTRE"}, \code{"TRE"}, \code{"GTRE_Z"} and \code{"TRE_Z"}, and, from 1.3.0, the closed-form likelihoods \code{"GTRE_FML"}, \code{"TFE_WMLE"}, \code{"FD"}, \code{"PL80"}, \code{"PL80_MVTN"}, \code{"BC92"}, \code{"K1990"} and \code{"K1990modified"}. This lets the log-likelihood be evaluated away from the optimum, which profiling requires; \code{\link{sfa_diagnostics}} uses it. \code{FALSE} by default, since the closure captures the data and the draws. From 1.3.0 the closure also accepts \code{per_obs = TRUE}, returning one log-likelihood contribution per \emph{firm} (firms are the independent units), so the score-based tools work on these fits: \code{\link{influence_sfa}}, \code{vcov(type = "bhhh")}, the \pkg{sandwich} methods, \code{\link{TIC}} and \code{\link{vuong}}. Their rows and sandwich standard errors are firm-level, i.e. clustered by firm. The other \code{model_name}s (\code{"TFE"}, \code{"SSFE"}, \code{"SSRE"}, \code{"SSCRE"}, \code{"CSS"}, \code{"LS"}, \code{"KSS"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"}) are not fit by maximum likelihood here, so they have no score matrix; for those the argument has no effect and \code{psfm()} warns.}
\item{kss_L}{\code{"KSS"} only. The dimension \eqn{L} of the factor space. \code{"auto"} (the default) selects it by the Bai and Ng (2002) \eqn{IC_{p2}} criterion; give an integer to fix it.}
\item{kss_smooth}{\code{"KSS"} only. The smoothing penalty applied to each firm's residual trajectory before the basis is extracted. \code{"auto"} (the default) chooses it by generalized cross-validation; \code{0} turns smoothing off, which reduces the model to an unrestricted factor model.}
\item{kss_L_max}{\code{"KSS"} only. Largest \eqn{L} the automatic criterion may choose. The criterion is also capped internally at \eqn{\lfloor T/2 \rfloor} (see \sQuote{Details}) and warns if it selects that cap. An explicit \code{kss_L} is not subject to the \eqn{\lfloor T/2 \rfloor} cap, only to \eqn{\min(T-1, n-1)}.}
\item{sim_type}{Sequence used to build the draws for the simulated-ML models listed under \code{halton_num}; ignored by every other model. \code{"halton"} (default) reproduces previous behaviour exactly. \code{"sobol"} accepts \code{sim_scrambling}; \code{"torus"} accepts \code{sim_prime}; \code{"uniform"} is pseudorandom, a baseline rather than a recommendation. The same options as \code{\link{sfm}}, built by the same code; see \sQuote{Details} there. Each firm receives its own block of the sequence, and \code{rand.gtre} randomizes it by a shift modulo 1 (by the seed, for \code{"sobol"} and \code{"uniform"}).}
\item{antithetics}{Logical. If \code{TRUE}, half of each firm's draws are taken from the sequence and the other half created as their mirror images. Default \code{FALSE}.}
\item{sim_burn}{Number of leading elements of the sequence discarded. Default \code{NULL}, meaning 1000.}
\item{sim_scrambling}{Integer 0--3, \code{sim_type = "sobol"} only, passed to \code{randtoolbox::sobol}.}
\item{sim_prime}{Prime for \code{sim_type = "torus"}. Default \code{NULL} uses the generator's own choice.}
\item{mundlak}{Mundlak adjustment terms: a formula naming the regressors whose \bold{within-firm means} should be added to the frontier, or \code{~ .} for every frontier regressor. \code{NULL} (the default) adds nothing. Only for models that treat the firm effect as random; \code{"TFE"}, \code{"TFE_WMLE"}, \code{"FD"} and \code{"SSFE"} refuse it, because they have already removed the effect. See \sQuote{Details}.}
}
\details{\code{psfm()} fits panel frontiers. Choose the estimator with \code{model_name}:

\tabular{lll}{
 \bold{\code{model_name}} \tab \bold{what it models} \tab \bold{use it when} \cr
 \code{"TRE"} \tab true random effects: firm heterogeneity plus transient inefficiency, by simulated ML \tab firm effects are not inefficiency \cr
 \code{"TRE_Z"} \tab \code{"TRE"} with determinants of \eqn{\sigma_u} \tab as \code{"TRE"}, and covariates drive inefficiency \cr
 \code{"GTRE"} \tab four components: heterogeneity, persistent and transient inefficiency; route set by \code{estimator} \tab both kinds of inefficiency matter \cr
 \code{"GTRE_Z"} \tab \code{"GTRE"} by simulated ML, with determinants \tab as \code{"GTRE"}, and covariates drive inefficiency \cr
 \code{"GTRE_FML"} \tab the four-component model by exact (closed-skew-normal) ML \tab a balanced panel; same as \code{"GTRE"}'s default \cr
 \code{"GTRE_SEQ1"} \tab four components in two steps, ML second stage \tab \eqn{T} is large; inconsistent for fixed \eqn{T} \cr
 \code{"GTRE_SEQ2"} \tab four components in two steps, moment second stage \tab \eqn{T} is large; inconsistent for fixed \eqn{T} \cr
 \code{"TFE"} \tab Greene's true fixed effects \tab firm effects may correlate with the regressors; long panels \cr
 \code{"TFE_WMLE"} \tab within ML true fixed effects (Chen, Schmidt and Wang 2014) \tab as \code{"TFE"}, on short panels \cr
 \code{"FD"} \tab first differences (Wang and Ho 2010) \tab fixed effects with inefficiency determinants \cr
 \code{"PL80"} \tab time-invariant inefficiency (Pitt and Lee 1980) \tab inefficiency is constant within a firm \cr
 \code{"BC92"} \tab inefficiency decaying exponentially (Battese and Coelli 1992) \tab inefficiency trends smoothly over time \cr
 \code{"PL80_MVTN"} \tab time-varying inefficiency, correlated within firm \tab a balanced panel; no imposed time path \cr
 \code{"K1990"} \tab logistic time path (Kumbhakar 1990) \tab a richer time path than \code{"BC92"} \cr
 \code{"K1990modified"} \tab quadratic time path \tab a richer time path than \code{"BC92"} \cr
 \code{"SSFE"} \tab Schmidt-Sickles within estimator; effect is inefficiency \tab no distribution for \eqn{u}; time-invariant \cr
 \code{"SSRE"} \tab Schmidt-Sickles GLS (random effects) \tab as \code{"SSFE"}, effects uncorrelated with regressors \cr
 \code{"SSCRE"} \tab \code{"SSRE"} with Mundlak means \tab as \code{"SSRE"}, effects may correlate \cr
 \code{"CSS"} \tab a firm-specific quadratic in time \tab no distribution; firms may cross over \cr
 \code{"LS"} \tab one common time pattern scaled per firm \tab no distribution; rankings fixed over time \cr
 \code{"KSS"} \tab a smooth time-varying effect with an estimated basis \tab no distribution; a balanced panel
}

The subsections below cover the panel index and formula, each model family, starting values, standard errors and efficiency prediction. The Monte Carlo measurements behind several of the caveats are collected in \code{vignette("sfa_numerical_notes", package = "sfa")}.

\subsection{Formula, panel index and time}{
Every model needs \code{individual}, the column that identifies panel units. The models with a time path also use \code{time}; without it, rows are taken to be in chronological order within each unit. A \code{plm::pdata.frame} supplies both.

For \code{"TRE_Z"} and \code{"GTRE_Z"} the segment after the first pipe holds the determinants of \eqn{\sigma_u}, and for \code{"GTRE_Z"} a second segment holds those of \eqn{\sigma_h}. These models place \eqn{z'\delta} on the \emph{variance}, \eqn{\sigma_u = \sqrt{\exp(z'\delta)}}, so \eqn{\delta} is twice what \code{\link{sfm}}'s default \code{z_link = "sd"} reports for the same data.

\code{mundlak} adds within-firm means to the frontier for the random-effects models; see \sQuote{Mundlak adjustment terms} below.
}

\subsection{The four-component model (GTRE)}{
\code{"GTRE"} and \code{"GTRE_FML"} fit
\deqn{y_{it} = \beta_0 + x_{it}'\beta + r_i - h_i + v_{it} - u_{it},}
with \eqn{r_i, v_{it}} two-sided and \eqn{h_i, u_{it}} half-normal, separating \emph{persistent} inefficiency \eqn{h_i} from \emph{transient} inefficiency \eqn{u_{it}} (Filippini and Greene 2016; Colombi, Kumbhakar, Martini and Vittadini 2014).

\bold{Choosing the route.} Since version 1.1.4 the route is chosen with \code{estimator}. \code{model_name = "GTRE"} now defaults to \code{estimator = "fiml"}; through 1.1.3 it fitted simulated ML, now \code{estimator = "sml"}. Scripts that pass \code{"GTRE"} therefore get a different estimator, and are warned once per call. \code{"GTRE_FML"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are unaffected.

\bold{Full information ML.} With \eqn{A = -[1_T, I_T]}, \eqn{V = \mathrm{diag}(\sigma_h^2, \sigma_u^2 I_T)} and \eqn{\Sigma = \sigma_v^2 I_T + \sigma_r^2 1_T 1_T'}, the firm-level density has a closed-skew-normal form, so nothing is simulated and repeated fits agree exactly. The \eqn{(T+1)}-dimensional normal CDF it needs has a diagonal-plus-rank-one covariance, so it reduces to a one-dimensional Gauss--Hermite integral. \code{"GTRE_FML"} reports raw standard deviations \code{sigr}, \code{sigv}, \code{sigh}, \code{sigu}, not the \eqn{\lambda}/\eqn{\sigma} parameterization of \code{"GTRE"} under \code{"sml"}. It requires a \bold{balanced} panel and errors otherwise.

\bold{Starting values.} The likelihood has a boundary optimum at \eqn{\sigma_h = 0}, where the model reduces to \code{"TRE"} and the intercept absorbs \eqn{E[h] = \sigma_h\sqrt{2/\pi}}. A fit can land there because it got stuck, or because that is the maximum for the sample (see the next subsection). The start search addresses only the first case. Following Colombi (2010) and Colombi, Martini and Vittadini (2011), \code{psfm()} evaluates two candidate starts, the random-effects decomposition and the \code{"GTRE_SEQ2"} estimate, and begins from the better one. \code{start_search} records the choice. \code{start_val} overrides it.

\bold{The persistent split is the weakest part of the model.} \eqn{\sigma_r} and \eqn{\sigma_h} enter the density only through \eqn{\sigma_r^2 + \sigma_h^2}, and their separation rests on the skewness term alone. Expect them to be the least precisely estimated parameters even at the global optimum.
}

\subsection{When a persistent scale collapses to zero}{
\code{"GTRE"} carries two persistent components, \eqn{\sigma_h} and \eqn{\sigma_r}, and a given sample cannot always separate them. When it cannot, one scale goes to zero and the other absorbs its variation. This is common, not exceptional, even when both components are present. \code{psfm} reports it in \code{$sigh_at_bound} and \code{$sigr_at_bound} and warns, naming the surviving scale.

\bold{A collapse is usually the right answer.} The boundary fit often has a higher likelihood than the true parameters, just as \eqn{\sigma_u = 0} is the maximum likelihood estimate under wrong skewness in a cross-section (Waldman 1982; see \code{\link{skewness_test}}). The start search above removes trapped optimizers; what remains is the estimate. Bounding either scale away from zero would corrupt exactly the samples where the boundary is the answer. Read the two persistent scales together and treat their split as unidentified in that sample; their sum is still informative.

\bold{What a collapse does to the other parameters.} It conserves the persistent variance and moves the persistent mean into the intercept, since \eqn{E[y] = \beta_0 - E[h] - E[u]}. The survivor does \emph{not} come back at \eqn{\sqrt{\sigma_h^2 + \sigma_r^2}}{sqrt(sigma_h^2 + sigma_r^2)}:

\tabular{lll}{
  \tab \bold{surviving scale} \tab \bold{intercept} \cr
  \eqn{\sigma_h \to 0}{sigma_h -> 0} \tab \eqn{\hat\sigma_r = \sqrt{\sigma_r^2 + \sigma_h^2(1-2/\pi)}}{sigma_r_hat = sqrt(sigma_r^2 + sigma_h^2 (1 - 2/pi))} \tab shifts by \eqn{-E[h]}{-E[h]} \cr
  \eqn{\sigma_r \to 0}{sigma_r -> 0} \tab \eqn{\hat\sigma_h = \sqrt{\sigma_h^2 + \sigma_r^2/(1-2/\pi)}}{sigma_h_hat = sqrt(sigma_h^2 + sigma_r^2 / (1 - 2/pi))} \tab shifts by \eqn{+\Delta E[h]}{+Delta E[h]}
}

So a collapsed fit displaces the intercept by a known, substantial amount, and any level read off \code{coef()} moves with it.

\bold{What a collapse says about the data.} A collapsed \eqn{\sigma_h} is wrong skewness between firms: average the OLS residuals within each firm and look at their skew. A positive between-firm skew predicts a \eqn{\sigma_h} collapse well. Treat this as a reading of the sign, not a test. The firm mean also contains \eqn{-\bar u_i}, which is itself skewed, so on short panels the between-firm skew leans negative even when \eqn{\sigma_h = 0}. \code{\link{skewness_test}} does not apply to the between-firm residuals as it stands.

\bold{Whether your panel can support the split.} The governing quantity is not \eqn{\sigma_h/\sigma_r} but the persistent variance of the firm-mean composite against the transient variance left in it:
\deqn{S/N = \frac{\sigma_r^2 + \sigma_h^2(1 - 2/\pi)}{[\sigma_u^2(1 - 2/\pi) + \sigma_v^2]/T}}{S/N = [sigma_r^2 + sigma_h^2 (1 - 2/pi)] / ([sigma_u^2 (1 - 2/pi) + sigma_v^2] / T)}
As a rule of thumb, below \eqn{S/N \approx 1}{S/N about 1} about half of fits merge the two components, and it takes roughly \eqn{S/N > 5} to get under 10\%. Two consequences:
\itemize{
  \item \bold{Raising \eqn{\sigma_h/\sigma_r} does not help}; it only decides which scale collapses. Shrinking \eqn{\sigma_r} lowers \eqn{S/N} and raises the collapse rate.
  \item \bold{More transient inefficiency hurts}: \eqn{\sigma_u} is in the denominator.
}

\bold{Firms or periods.} \eqn{S/N} covers \eqn{T} but not \eqn{N}, and the two help through different channels. Under the default \code{estimator = "fiml"}, adding firms reduces collapse more than adding periods. Under \code{"sml"} the reverse holds, and a short panel stays a dead end however many firms it has. So if an unbalanced panel forces \code{"sml"}, panel length is what buys the split. Which estimator collapses less often depends on the design: where the random effect dominates, \code{"sml"} can collapse less often than \code{"fiml"}.
}

\subsection{Two-step estimators (GTRE_SEQ1, GTRE_SEQ2)}{
Both take the random-effects decomposition from \code{plm} and pass its two
output series to a second stage: \code{"GTRE_SEQ1"} fits an intercept-only
normal--half-normal likelihood to each, \code{"GTRE_SEQ2"} inverts their second
and third moments. The second stage treats those series as draws from
\eqn{v_{it}-u_{it}} and \eqn{r_i-h_i}, which at finite \eqn{T} they are not:
\code{alpha_hat} is the shrunken predictor \eqn{c\,(\alpha_i+\bar e_i)} with
\eqn{c = T\sigma_\alpha^2/(T\sigma_\alpha^2+\sigma_\epsilon^2)}, and
\code{eps_hat} is the quasi-demeaned residual
\eqn{(1-\theta)\alpha_i + e_{it} - \theta\bar e_i} with
\eqn{\theta = 1-\sigma_\epsilon/\sqrt{T\sigma_\alpha^2+\sigma_\epsilon^2}}.

This attenuates the third central moment on the \eqn{\epsilon} side by
\eqn{(1-\theta/T)^3+(T-1)(-\theta/T)^3} (about 0.81 at \eqn{T=10}) and scales
the \eqn{\alpha} side by \eqn{c^2} and \eqn{c^3}. None of these factors tends to
1 as \eqn{N} grows with \eqn{T} fixed, so \bold{both estimators are inconsistent
for fixed \eqn{T}}, and each call warns. Raising \eqn{T} removes the bias;
raising \eqn{N} does not. For a consistent estimator use
\code{model_name = "GTRE"}.

Their delta-method standard errors use \emph{iid} moment variances, so they
describe variability around each estimator's probability limit, not around the
truth. They are also optimistic on the \eqn{\epsilon} side, whose residuals are
correlated within firm.
}

\subsection{Simulated maximum likelihood (TRE, TRE_Z, GTRE_Z, GTRE with sml)}{
These models integrate the firm effects by simulation over \code{halton_num} draws per firm, from the sequence \code{sim_type} selects (Filippini and Greene 2016). \code{"TRE_Z"} and \code{"GTRE_Z"} add determinants of the one-sided scales. The simulated log-likelihood depends on the draws: \code{rand.gtre} changes the objective itself, so log-likelihoods from different draws are not comparable. Use \code{\link{simulation_se}} to measure simulation error and raise \code{halton_num} to reduce it. \code{"TRE_Z"} becomes unstable on small panels (below about 60 firms).
}

\subsection{True fixed effects and first differences}{
\code{"TFE"} and \code{"TFE_WMLE"} fit the same model,
\deqn{y_{it} = \alpha_i + x_{it}'\beta + v_{it} - u_{it},\quad v \sim N(0,\sigma_v^2),\ u \sim N^+(0,\sigma_u^2),}
with \eqn{\alpha_i} an unrestricted individual effect. They differ in how they treat \eqn{\alpha_i}.

\code{"TFE"} is Greene's (2005) true fixed effects estimator: the composed-error likelihood with the \eqn{N} intercepts estimated jointly. The package concentrates each \eqn{\alpha_i} out at every evaluation (the unique root of a strictly decreasing score, solved by safeguarded Newton), so the parameter vector and \code{out} layout match \code{"TFE_WMLE"}. It suffers the \emph{incidental parameters problem}: with \eqn{T} fixed, \eqn{\sigma_u} is biased upward and efficiency scores are too pessimistic.

\code{"TFE_WMLE"} is Chen, Schmidt and Wang's (2014) within maximum-likelihood estimator, which applies ML to deviations from individual means. The \eqn{\alpha_i} drop out, so it has no incidental parameters problem. The firm effects are recovered afterwards from the mean residual plus \eqn{\sqrt{2/\pi}\,\sigma_u}.

Caveats:
\itemize{
  \item \bold{The two log-likelihoods are not comparable}: \code{"TFE"} is the likelihood of the data, \code{"TFE_WMLE"} of the within deviations. Do not use \code{logLik()}, \code{AIC()} or \code{BIC()} to choose between them.
  \item \bold{\code{"TFE"} changed meaning.} Through sfa 1.1.2 it selected the Chen-Schmidt-Wang estimator, now \code{"TFE_WMLE"}. Scripts that pass \code{"TFE"} get Greene's estimator and are warned once per call.
  \item \bold{Greene's likelihood is degenerate.} Its supremum is always on the \eqn{\sigma_v \to 0} boundary: setting \eqn{\alpha_i = \max_t (y_{it} - x_{it}'\beta)} makes every composed error weakly negative, and the deterministic-frontier likelihood that remains is finite. An unconstrained search that leaves the interior basin returns \eqn{\sigma_v \approx 0} with meaningless standard errors, and with large \eqn{\lambda} and small \eqn{T} no interior maximum need exist. \code{tfe_lambda_max} bounds the search, and a fit at the bound warns.
  \item \bold{Both problems shrink as \eqn{T} grows.} On short panels \code{"TFE"} pins at the bound; on long ones it converges and \code{tfe_lambda_max} is irrelevant. \code{"TFE_WMLE"} is unaffected by \eqn{T} either way.
}

\code{"FD"} is the first-difference estimator of Wang and Ho (2010), which also removes the fixed effect and lets determinants scale inefficiency.
}

\subsection{Error-components frontiers (PL80, BC92, K1990, PL80_MVTN)}{
\code{"PL80"} (Pitt and Lee 1980) and \code{"BC92"} (Battese and Coelli 1992) fit \eqn{y_{it} = x_{it}'\beta + v_{it} - B_{it}u_i} (production; the sign flips for cost), with \eqn{v_{it}} iid normal and \eqn{u_i} half-normal and constant within a firm. \eqn{B_{it} = 1} for \code{"PL80"} and \eqn{B_{it} = \exp(-\eta(t - T_{ref}))} for \code{"BC92"}, with \eqn{T_{ref}} the last period in the panel. The closed-form likelihood integrates \eqn{u_i} out of the joint density, and nests the cross-sectional normal-half-normal likelihood when every firm has one observation. \code{time} supplies each observation's period.

\code{"K1990"} uses \eqn{B_{it} = (1 + \exp(b t + c t^2))^{-1}} and \code{"K1990modified"} uses \eqn{B_{it} = 1 + d(t - T_i) + e(t - T_i)^2}. All four share one likelihood, differing only in \eqn{B_{it}}. \code{"K1990"}'s \code{b} and \code{c} are weakly identified: quite different pairs give near-identical paths, so interpret the fitted \eqn{B_{it}} path rather than either coefficient.

\code{"PL80_MVTN"} drops the common-\eqn{u_i} restriction and lets inefficiency vary over time with within-firm correlation; see \code{\link{PL80_MVTN}}.
}

\subsection{Schmidt-Sickles and time-varying fixed effects}{
These estimators read inefficiency off a firm effect, as distance from the best firm, and assume no distribution for \eqn{u}. None returns a log-likelihood, so \code{logLik()}, \code{AIC()} and \code{BIC()} are not defined for them (see \code{\link[sfa]{sfareg-methods}}).

\bold{\code{"SSFE"}, \code{"SSRE"}, \code{"SSCRE"}.} The within, GLS and correlated-random-effects members of the Schmidt and Sickles (1984) family. \code{"SSFE"} is the LSDV estimator; it attributes all within-firm deviation to inefficiency, unlike \code{"TFE"}, which keeps a separate noise term. \code{"SSRE"} is more efficient and identifies coefficients on time-invariant regressors, but requires the effects to be uncorrelated with the regressors. In a frontier that says a firm's inefficiency is unrelated to its input choices, which is rarely credible.

\code{"SSCRE"} models the correlation by adding the within-firm means of the time-varying regressors (Mundlak 1978). The slopes on the original regressors then equal the within estimates exactly, and time-invariant regressors stay identified. The means' coefficients carry a \code{.mean_} prefix; testing them jointly against zero is the Mundlak form of the Hausman test. A regressor constant within firm gets no mean term. The variance components (and so \eqn{\theta} and the BLUPs) are the Swamy--Arora components of the model \emph{without} the means. On an unbalanced panel \pkg{plm}'s Swamy--Arora step is singular for the augmented model, so \code{"SSCRE"} obtains the same components through \pkg{plm}'s Wallace--Hussain step and checks that they agree.

\bold{\code{"CSS"} and \code{"LS"}.} Both let the firm effect vary over time. Cornwell, Schmidt and Sickles (1990) give each firm its own quadratic, \eqn{\alpha_{it} = \theta_{i1} + \theta_{i2}t + \theta_{i3}t^2}: \eqn{3N} parameters, and firms may overtake one another. Lee and Schmidt (1993) impose one common pattern scaled per firm, \eqn{\alpha_{it} = \delta_t\alpha_i}: \eqn{N + T - 1} parameters, and the ranking of firms cannot change. Crossover or no crossover is the substantive choice between them.

\code{"CSS"} residualizes \eqn{y} and \eqn{X} on each firm's \eqn{(1, t, t^2)}. A firm needs \eqn{T_i \ge 4} to contribute to \eqn{\beta}; with three periods its quadratic fits exactly and \code{sfa} warns. A regressor spanned by those quadratics, such as a time trend, is refused as a rank failure. \code{"LS"} is a rank-one factor model fitted by alternating least squares, which handles unbalanced panels; \eqn{\delta_1 = 1} fixes the scale, with a norm constraint as fallback when the first period's effect is near zero.

For both, inefficiency is measured against the best firm \emph{in each period}, \eqn{u_{it} = \max_j \alpha_{jt} - \alpha_{it}}. The intercept goes into the firm effect, so \code{coef()} returns slopes only.
}

\subsection{Kneip, Sickles and Song (KSS)}{
\code{"KSS"} writes the firm effect as \eqn{\alpha_{it} = \sum_{r=1}^{L}\theta_{ir} g_r(t)}, estimating both the basis functions \eqn{g_r} and the dimension \eqn{L}. The other time-varying estimators are special cases:

\tabular{lll}{
  \strong{model} \tab \strong{\eqn{L}} \tab \strong{basis} \cr
  \code{"SSFE"} \tab 1 \tab fixed, constant in \eqn{t} \cr
  \code{"LS"}   \tab 1 \tab free \cr
  \code{"CSS"}  \tab 3 \tab fixed to \eqn{\{1, t, t^2\}} \cr
  \code{"KSS"}  \tab estimated \tab estimated
}

Estimation follows the original three steps, iterated with \eqn{\beta} until the residual sum of squares settles: centre by period (which imposes \eqn{\sum_i u_i(t) = 0}), smooth each firm's residual trajectory and take the eigendecomposition of their covariance to get \eqn{\hat g_r}, then fit loadings by least squares. Smoothing is what makes the basis a set of smooth functions rather than arbitrary period effects.

Caveats:
\itemize{
  \item \bold{\eqn{L} is chosen by Bai and Ng's (2002) \eqn{IC_{p2}}, not the threshold rule of the original paper}, and a different rule can choose a different \eqn{L}. Set \code{kss_L} to fix it, and inspect \code{$kss$eigenvalues} to see whether the choice was close. \code{kss_smooth} defaults to generalized cross-validation.
  \item \bold{The criterion needs a long enough panel.} On a short panel the factors can span almost the whole time dimension and no penalty competes, so an uncapped criterion picks its maximum. The automatic search is therefore capped at \eqn{\lfloor T/2 \rfloor}, which recovers the true rank from about \eqn{T = 8}. When the choice lands on the cap, \code{psfm} warns: read it as \dQuote{at least this many}, and prefer a longer panel to a larger \code{kss_L_max}.
  \item \bold{\code{alpha_hat} is not comparable across models.} Period centering removes the period mean from \code{"KSS"}'s \code{alpha_hat} but not from \code{"CSS"}'s or \code{"LS"}'s. Compare \code{u_hat}, a within-period contrast that is invariant to the shift.
}
}

\subsection{Mundlak adjustment terms}{
The true random effects model assumes the firm effect is independent of the regressors. In production data it usually is not, and the effect then leaks into the slopes and inefficiency. Mundlak's (1978) device models the correlation, \eqn{\alpha_i = \pi'\bar{x}_i + \delta_i}. Substituting into the frontier leaves the same likelihood with a wider design matrix, so this is an argument, not a new \code{model_name}.

Karagiannis and Kellermann (2019) extend it to \eqn{\alpha_i = \pi'\bar{x}^*_i + \gamma'z_i + \delta_i}, with environmental variables \eqn{z_i} beyond the producer's control. Those enter the formula as ordinary time-invariant regressors. \code{mundlak} builds the within-firm means against the same panel index the fit uses. Naming a \emph{subset} of the regressors gives the Hausman-Taylor variant.

\code{mundlak} removes the bias when the effect is correlated with the regressors, but costs precision when it is not. It is not a default.

A variable constant within firm gets no mean column, since it would equal its own mean. Every \code{mundlak} fit warns of collinearity in the \emph{starting-value} regression; this is expected, because \eqn{x} and \eqn{\bar{x}} coincide in a between-individual regression. The starts fall back to pooled OLS and the full formula is estimated.
}

\subsection{Starting values and optimization}{
The maximum-likelihood models run a staged optimizer (nlminb or bobyqa, optionally \code{psoptim}, then \code{optim} with \code{"L-BFGS-B"}). The defaults (\code{maxit.bobyqa}/\code{maxit.psoptim}/\code{maxit.optim} of 5000/100/1000) were raised in 1.1.3 from 100/10/10, which stopped several models short of the optimum. They stay below \code{\link{sfm}}'s because the simulated-ML models are expensive to evaluate. Lower them for a quick approximate fit.

Random-effects starting values come from a between regression; see \code{collinear_action} when it is rank deficient. \code{"GTRE_FML"} runs its own start search (see above).

\bold{Profiling.} With \code{keep_objective = TRUE} the fit carries \code{$objective}, the negative log-likelihood the optimizer minimized. Evaluating it away from the optimum separates a weakly identified parameter from a badly estimated one:

\preformatted{fit <- psfm(y ~ x1 + x2, "GTRE", d, individual = "id",
            estimator = "sml", keep_objective = TRUE)
ih   <- which(names(coef(fit)) == "sigh")
grid <- seq(0.05, 1, by = 0.05)
ll   <- vapply(grid, function(v) \{
          p <- fit$opt$par; p[ih] <- v; -fit$objective(p) \}, numeric(1))
plot(grid, ll - max(ll), type = "l")}

A flat curve means no estimator could recover the parameter from these data; a sharp one means the difficulty lies elsewhere.
}

\subsection{Standard errors}{
Standard errors come from the numerical Hessian and need \code{optHessian = TRUE}. \code{OPG_calc = TRUE} adds outer-product-of-gradients errors (\code{out_opg}, \code{out_sandwich} for \code{"GTRE_Z"}). With \code{keep_objective = TRUE}, \code{vcov(type = "bhhh")} and the \pkg{sandwich} methods give firm-clustered errors for the maximum-likelihood models. For simulated-ML fits, \code{\link{simulation_se}} measures the part of the uncertainty due to the draws, and \code{\link{psfm_bootstrap}} gives parametric bootstrap errors for the GTRE family. The two-step estimators' delta-method errors are discussed above. At a boundary fit the affected standard errors are not meaningful.
}

\subsection{Efficiency prediction}{
The GTRE family returns transient efficiency \code{U} per observation and persistent efficiency \code{H} per firm. \code{"TRE"}-type models return \code{U}. \code{"TFE"}, \code{"TFE_WMLE"}, \code{"FD"}, \code{"PL80"} and \code{"BC92"} return \code{exp_u_hat}, \eqn{E[\exp(-u)|\epsilon]} (Battese and Coelli 1988), generalized to the panel case for the error-components models by pooling a firm's observations and scaling by \eqn{B_{it}}. The distribution-free estimators return \code{u_hat} and \code{exp_u_hat} measured from the best firm. \sQuote{Value} gives the details per model; \code{\link{efficiency}} gives a common interface.
}
}
\references{Filippini and Greene (2016, JPA); Colombi, Kumbhakar, Martini and Vittadini (2014, JPA) on the closed-skew-normal four-component likelihood, with Colombi (2010) and Colombi, Martini and Vittadini (2011) on the two-step estimator used to seed it; Greene (2005, JoE and JPA); Wang and Ho (2010, JoE); Chen, Schmidt and Wang (2014, JoE); Kumbhakar (1990, JoE); Pitt and Lee (1980, JoE); Battese and Coelli (1992, JPA); Battese and Coelli (1988, JoE); Jondrow, Lovell, Materov and Schmidt (1982, JoE); Schmidt and Sickles (1984, JBES)}
\author{David H. Bernstein}

\seealso{\code{\link{sfm}} for cross-sectional models, \code{\link{psfm_bootstrap}} for bootstrap standard errors, \code{\link{data_gen_p}} for simulating panel data with known true parameters, \code{\link{sfareg-methods}} for extracting results, and \code{vignette("sfa_numerical_notes", package = "sfa")} for the measurements behind the caveats above.}

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
\item{inefdec}{The orientation this fit used, as a scalar logical: \code{TRUE} for a production frontier, \code{FALSE} for a cost one. Recorded at fit time so that post-estimation functions (\code{\link{efficiency}} with \code{logDepVar = FALSE}, the robust accessors, the composed CDF) do not have to re-evaluate the \code{inefdec} argument of the recorded call, which fails when it was passed as a variable that is no longer in scope.}
\item{data}{The data used in estimation.}

\bold{Additional model-specific components:}

\bold{For GTRE_FML model:}
\item{start_search}{A list recording the starting-value search: \code{n_tried}, the log-likelihood at each candidate (\code{loglik}), and which was \code{chosen} (\code{"random-effects"} or \code{"two-step"}). \code{NULL} when \code{start_val} was supplied.}
\item{stationarity}{List, present when \code{optHessian = TRUE} and none of the four scales is on its lower bound (a collapsed scale is reported by \code{sigh_at_bound}/\code{sigr_at_bound} instead): the Newton-decrement check described under \code{stationarity} in \code{\link{sfm}}. When the decrement exceeds \code{1e-4}, the fit is polished with \code{nlminb()}, the four scales on the log scale, and the result is kept only if the log-likelihood rises; \code{polished}, \code{gain} and \code{decrement_after} record what happened.}
\item{U}{Predicted transient (time-varying) technical efficiency for each observation.}
\item{H}{Predicted persistent (time-invariant) technical efficiency for each individual, named by individual.}

\bold{For GTRE and GTRE_Z models:}
\item{H}{Predicted time-invariant technical efficiency for each individual.}

\bold{For GTRE_Z model:}
\item{out_opg}{OPG (outer product of gradients) standard errors and t-values along with parameters.}
\item{out_sandwich}{sandwich standard errors and t-values along with parameters.}

\bold{For GTRE, GTRE_Z, TRE and TRE_Z models:}
\item{U}{Predicted time-varying technical efficiency for each observation.}

\bold{For TFE and TFE_WMLE models:}
\item{r_hat_m}{Estimated individual-specific fixed effects, named by individual. For \code{"TFE"} these are the maximum-likelihood estimates of \eqn{\alpha_i} directly; for \code{"TFE_WMLE"} they are recovered after estimation from the mean residual plus \eqn{\sqrt{2/\pi}\,\sigma_u}, since the within likelihood does not identify them.}
\item{exp_u_hat}{Predicted technical efficiency, \eqn{E[\exp(-u)|\epsilon]} (Battese and Coelli, 1988).}

\bold{For TFE model only:}
\item{u_hat}{Predicted inefficiency, \eqn{E[u|\epsilon]} (Jondrow et al., 1982).}

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

\bold{For the GTRE model:}
\item{sigh_at_bound, sigr_at_bound}{\code{TRUE} when that persistent scale has collapsed to the zero boundary, judged relative to \code{sigma}. Often the correct maximum likelihood estimate rather than a failure; see \sQuote{Details} -- so it is reported rather than suppressed, and \code{psfm} warns naming the other persistent scale, since the two must be read together.}

\bold{For CSS and LS models:}
\item{alpha_hat}{The fitted firm effect for every observation, \eqn{\alpha_{it}} -- time-varying, unlike \code{SSFE}'s one value per firm.}
\item{u_hat, exp_u_hat}{Inefficiency and efficiency per observation, measured against the best firm in the same period.}
\item{frontier_alpha}{The best firm effect in each period, i.e. the benchmark \code{u_hat} is measured from.}
\item{sigma_v, residuals, df.residual}{The residual standard deviation, the residuals net of the firm effect, and the residual degrees of freedom (\eqn{NT - K - 3N} for \code{CSS}, \eqn{NT - K - (N + T - 1)} for \code{LS}).}
\item{theta}{\code{CSS} only. An \eqn{N \times 3} matrix of the firm-specific quadratic coefficients.}
\item{delta, alpha_i, ls_iterations, ls_converged}{\code{LS} only. The common temporal pattern (with \code{delta[1] == 1}), the firm-specific scales, and the state of the alternating least-squares iteration.}

\bold{For SSRE and SSCRE models:}
\item{alpha_hat}{The BLUPs of the firm effects, one per firm. Under \code{SSCRE} these are net of the part the Mundlak means explain.}
\item{theta, ercomp}{The GLS transformation parameter and the full variance-components decomposition from \code{\link[plm]{ercomp}}.}
\item{mundlak_terms}{\code{SSCRE} only. The names of the added group-mean regressors.}
\item{plm_fit}{The underlying \code{\link[plm]{plm}} object, for users who want its own methods.}

\bold{For the KSS model:}
\item{kss}{A list: the selected dimension \code{L}, the smoothing penalty \code{kappa}, the estimated basis (a \eqn{T \times L} matrix of orthonormal columns), the \eqn{n \times L} matrix of \code{loadings}, the full \code{eigenvalues} sequence the dimension was chosen from, and the iteration state.}

\bold{For PL80 and BC92 models:}
\item{exp_u_hat}{Predicted technical efficiency for each observation, E[exp(-B_it*u_i)|epsilon_i] (generalizes Battese and Coelli's (1988, JoE) point estimator to the panel case, pooling all of a firm's observations to predict u_i and scaling by B_it per period).}
}

\note{
Standard errors require \code{optHessian = TRUE}.

The GTRE_SEQ1 and GTRE_SEQ2 models use sequential estimation methods and do not return optimization objects or starting values. The SSFE, CSS, LS, SSRE, SSCRE and KSS models are likewise not maximum-likelihood based and return no optimization object, starting values, or log-likelihood. All panel models require the \code{individual} argument to identify panel units.
}
\examples{
\donttest{
library(sfa)     

## A deliberately small panel with few simulation draws, so the examples
## build quickly. Do not read the estimates as a serious fit. TRE_Z becomes
## unstable below about 60 firms, so this is close to the floor.
data_trial <- data_gen_p(t=6, N=60, rand = 100,
                         sig_u = 1,  sig_v = 0.3,
                         sig_r = .2, sig_h = .4,
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

## rand.psoptim is not optional with PSopt = TRUE: the particle swarm draws
## from the session RNG, so results otherwise change between builds.
max_tre_z   <-  psfm(formula    = y_tre_z ~ x1 +x2| z_gtre,
                     model_name = "TRE_Z",
                     data       = data_trial,
                     individual = "name",
                     halton_num = 50,
                     rand.psoptim = 7L,
                     PSopt      = TRUE)

## Time-varying fixed effects, no distribution assumed for u. CSS gives each
## firm its own quadratic in time; LS gives all firms one common pattern.
css_fit     <-  psfm(formula = y_tre ~ x1 + x2, model_name = "CSS",
                     data = data_trial, individual = "name", time = "year")
ls_fit      <-  psfm(formula = y_tre ~ x1 + x2, model_name = "LS",
                     data = data_trial, individual = "name", time = "year")
ls_fit$delta      # the common temporal pattern, normalized to delta[1] = 1

## The Schmidt-Sickles family: within, GLS, and correlated random effects.
## SSCRE's slopes equal SSFE's within slopes exactly (Mundlak 1978); the
## .mean_ coefficients are the Mundlak form of the Hausman test.
ssre_fit    <-  psfm(formula = y_tre ~ x1 + x2, model_name = "SSRE",
                     data = data_trial, individual = "name", time = "year")
sscre_fit   <-  psfm(formula = y_tre ~ x1 + x2, model_name = "SSCRE",
                     data = data_trial, individual = "name", time = "year")

## Kneip, Sickles and Song (2012): the temporal basis and its dimension are
## both estimated. Inspect what the criterion chose before trusting it.
kss_fit     <-  psfm(formula = y_tre ~ x1 + x2, model_name = "KSS",
                     data = data_trial, individual = "name", time = "year")
kss_fit$kss$L                 # selected dimension
kss_fit$kss$eigenvalues       # was the choice close?
}
}
