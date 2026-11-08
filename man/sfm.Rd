\name{sfm}
\alias{sfm}
\title{Cross-Sectional Stochastic Frontier Models}
\description{Implementation of the cross-sectional stochastic frontier model across an array of distributional assumptions for both v and u (user specified). For panel models, see the psfm() call.}
\usage{sfm(formula, model_name = c("NHN", "NHN_Z", "NE", "NE_Z", "NR",
    "THT", "NTN", "NG", "NNAK", "NU", "NGE", "NLN", "NW", "tHN",
    "TSL", "NB", "NGB2"), data, maxit.bobyqa = 10000, maxit.nlminb = 500,
    maxit.psoptim = 1000, maxit.optim = 1000, REPORT = 1, trace = 2,
    pgtol = 0, start_val = FALSE, PSopt = FALSE, use.nlminb = "auto",
    use.bobyqa = "auto", optHessian = TRUE, inefdec = TRUE, upper =
    NA, Method = "L-BFGS-B", robust = c("mle", "mlqe", "psi",
    "mdpd"), c_mlqe = 0.2, eta = 0.01, alpha = 0.2, verbose = FALSE,
    start_from = NULL, Nsim = "auto", z_link = c("sd", "var"), vhet =
    NULL, uhet = NULL, muhet = NULL, scaling = NULL, shapehet = NULL,
    weights = NULL, wscale = TRUE, sim_type = c("halton", "sobol",
    "torus", "uniform"), antithetics = FALSE, sim_burn = NULL,
    sim_scrambling = 0L, sim_prime = NULL, sim_seed = NULL,
    rand.psoptim = NULL, keep_objective = FALSE, estimator = c("mle",
    "cols", "mols", "acols", "cmle"), cols_boot = 0, rand.cols = NULL,
    integration = c("quadrature", "simulate"))
}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{Which distributional pair to fit. The noise term \eqn{v} is
normal except where stated, and the one-sided inefficiency term \eqn{u} follows:
\describe{
  \item{\code{"NHN"}}{half-normal.}
  \item{\code{"NHN_Z"}}{half-normal with heteroskedastic \eqn{\sigma_u}, parameterized by the determinants \eqn{z} given in the second part of the formula (\code{y ~ x | z}) as \eqn{\sigma_u = \exp(z'\delta)}.}
  \item{\code{"NE"}}{exponential.}
  \item{\code{"NE_Z"}}{exponential with heteroskedastic \eqn{\sigma_u}, parameterized as for \code{"NHN_Z"}.}
  \item{\code{"NTN"}}{truncated normal.}
  \item{\code{"NR"}}{Rayleigh. A distinct family from \code{"NHN"}, not a reparameterization of it; see \sQuote{Details}.}
  \item{\code{"NU"}}{uniform.}
  \item{\code{"NG"}}{gamma. The likelihood determines \eqn{E[u] = \mu\sigma_u} far better than either factor; see \code{ng_starts} under \sQuote{Value}.}
  \item{\code{"NNAK"}}{Nakagami, with shape \eqn{m}. Nests \code{"NR"} at \eqn{m = 1}. The shape is weakly identified, so \eqn{\sigma_u} and \eqn{m} are strongly correlated across samples.}
  \item{\code{"NGE"}}{generalized exponential.}
  \item{\code{"NLN"}}{lognormal. The composed density is computed by deterministic quadrature; see \code{integration} and \sQuote{Details}.}
  \item{\code{"NW"}}{Weibull. Computed the same way as \code{"NLN"}.}
  \item{\code{"tHN"}}{half-normal, with Student-\emph{t} noise.}
  \item{\code{"TSL"}}{truncated skew-Laplace.}
  \item{\code{"THT"}}{half \emph{t}, with Student-\emph{t} noise.}
  \item{\code{"NGB2"}}{generalized beta of the second kind (Makiela and Mazur 2022), with scale \eqn{\sigma_u} and three shape parameters \eqn{(\nu, \psi, \tau)}. It nests \code{"NHN"}, \code{"NE"}, \code{"NG"} and \code{"NW"} as limits and the half-Student \eqn{t} exactly. Evaluated by deterministic quadrature; see \sQuote{Details}.}
  \item{\code{"NB"}}{binomial (Carree 2002). \bold{Corrected OLS only}: call it with \code{estimator = "cols"}; \code{sfm()} refuses it under maximum likelihood. The only inefficiency distribution here that can be skewed \emph{either} way; see \sQuote{Details}.}
} Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT, trace, pgtol, upper}{Accepted but not used by \code{sfm()}: the optimizer's own tracing is controlled by \code{verbose}, and its bounds are set internally.}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T or F)}
\item{maxit.nlminb}{Maximum number of iterations for the nlminb optimization routine.}
\item{use.nlminb}{\code{"auto"} (default), or a logical. Whether to run the \code{stats::nlminb} stage, which precedes bobyqa and uses an analytic gradient where one exists. Under \code{"auto"} it runs for \code{NHN}, \code{NE}, \code{NTN} and \code{NU}, where it reaches an optimum at least as good as bobyqa's, faster. \code{NHN_Z} and \code{NE_Z} also take nlminb under \code{"auto"}, starting from an unconstrained fit of the same model with \eqn{\log\sigma_v} free (the \code{vhet} path with an intercept only): an uncentred \eqn{z} puts a ridge between the intercept and slope of \eqn{\delta}, and bobyqa does not climb it within its evaluation budget.}
\item{use.bobyqa}{\code{"auto"} (default), or a logical. Whether to run the derivative-free bobyqa stage. Under \code{"auto"} it runs for every model \emph{except} those handled by nlminb above. For \code{NR} and \code{NGE} the nlminb path reaches a worse optimum, so bobyqa stays on for them.}
\item{optHessian}{Logical. Finish with \code{optim()} and a numerically differentiated Hessian? Standard errors and t-values require \code{TRUE}, the default. With \code{FALSE} the final \code{optim} stage is not run at all, so the estimates reported are the preceding stage's and can differ from the default fit's wherever that stage had not already converged.}
\item{inefdec}{Production or cost function}
\item{Method}{Method passed to \code{optim} for the final stage. Must be
\code{"L-BFGS-B"}: it is the only \code{optim} method that applies the lower
and upper bounds this package uses to keep scale parameters positive. Any
other value is an error.}
\item{robust}{One of \code{"mle"} (default), \code{"mlqe"}, \code{"psi"}, or
\code{"mdpd"}. Selects a robust divergence-based objective in place of
ordinary maximum likelihood; see \sQuote{Details}. Implemented for
\code{model_name = "NHN"} only; any other model errors if a
non-\code{"mle"} value is requested.}
\item{c_mlqe}{Numeric, default \code{0.20}. The MLqE power parameter
(\code{c} in \code{f^c}; the corresponding \code{q = 1-c}), used only when
\code{robust = "mlqe"}.}
\item{eta}{Numeric, default \code{0.01}. The Psi-likelihood power parameter,
used only when \code{robust = "psi"}.}
\item{alpha}{Numeric, default \code{0.2}. The MDPD (minimum density power
divergence) power parameter, used only when \code{robust = "mdpd"}.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{weights}{Optional observation weights, one per row of the data actually fitted (after missing-data handling; a length mismatch is an error, not recycled). They multiply the per-observation log-likelihood contributions. For \bold{frequency} weights this is exact: a weight of 2 gives the same estimates as including the row twice. For \bold{sampling} weights it gives a pseudo-likelihood whose point estimates are consistent but whose Hessian-based standard errors are not; use the \pkg{sandwich} methods this package registers, e.g. \code{sandwich::sandwich(fit)}. Cannot be combined with the robust divergence estimators, which reweight observations themselves.}
\item{wscale}{Rescale \code{weights} to sum to the number of observations. This puts the log-likelihood on the same scale as an unweighted fit, so \code{AIC}/\code{BIC} stay comparable; it does not move the estimates. Set \code{FALSE} for genuine frequency weights, where the unscaled total is the quantity of interest.}
\item{shapehet}{\bold{Experimental.} Optionally, a one-sided formula parameterizing the \emph{shape} of \eqn{u} rather than its scale, as \eqn{\mu_i = \mu\exp(z_i'\delta)}. Implemented for \code{model_name = "NG"} and \code{"NNAK"}, the two families with a shape parameter. It raises the log-likelihood sharply on data with a true shape effect, \bold{but its coefficients do not recover well in testing}, and \eqn{\sigma_V} tends to collapse to its lower bound, the signature of a degenerate optimum. Plain \code{"NG"} does not have this problem. Treat the coefficients as exploratory, check \eqn{\sigma_V} for a boundary collapse before reading them, and prefer \code{uhet} when you mean scale heterogeneity.}
\item{scaling}{Optionally, a one-sided formula of covariates imposing the \bold{scaling property} of Wang and Schmidt (2002) and Alvarez, Amsler, Orea and Schmidt (2006): \eqn{u_i = h(z_i,\delta)\,u^*_i} with \eqn{h = \exp(z_i'\delta)}, so the covariates scale a common draw \eqn{u^*}. Both \eqn{\sigma_U} and the pre-truncation mean move by the \emph{same} factor, which fixes the shape of the distribution and frees only its scale. Implemented for \code{model_name = "NTN"} only: for a half-normal \eqn{u}, \eqn{h(z)|N(0,\sigma^2)|} \emph{is} \eqn{|N(0,(h\sigma)^2)|}, so \code{uhet} already fits that model. Cannot be combined with \code{uhet} or \code{muhet}, which would undo the restriction. The reported \code{scale.} coefficients are \eqn{\delta}; \code{Zu.(Intercept)} and \code{Zmu.(Intercept)} are the scale and pre-truncation mean of the common draw.}
\item{start_from}{Optionally, a previously fitted \code{"sfareg"} object to take starting values from. Parameters are matched \bold{by name}, not by position: \code{"NE"} reports \code{sigv}/\code{sigu} where \code{"NHN"} reports \code{lambda}/\code{sigma}, so a positional copy would put a variance ratio into a standard deviation. Unmatched parameters keep the package's own start. Fitting a simpler model first and cascading is the remedy for awkward likelihoods such as \code{"NG"} (seeded from \code{"NE"}). The seeded fit is rarely better, but when it is the gain can be large, and nothing in the default fit warns that it stopped at a worse optimum.}
\item{Nsim}{Total number of simulation draws per observation for \code{"NLN"} and \code{"NW"} under \code{integration = "simulate"}; ignored otherwise, with a warning if set under quadrature. \code{"auto"} (the default) uses \code{max(200, ceiling(3*sqrt(n)))}. Half the draws come from each of the two proposals described under \sQuote{Details}. Simulated ML is consistent only if the draw count grows with the sample size: at a fixed \code{Nsim} every parameter except the frontier slopes converges to the wrong value. A numeric value below the automatic floor warns. \code{"NLN"} and \code{"NW"} need similar draw counts.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
  \item{keep_objective}{Logical. Store the likelihood function on the fitted object, so that \code{\link{sfa_diagnostics}} can compute the gradient at the optimum and draw likelihood slices. Defaults to \code{FALSE}: the stored closure carries the estimation data, so a saved fit can be many times larger. Everything else \code{sfa_diagnostics()} reports works without it.}
\item{estimator}{\code{"mle"} (default), or \code{"cols"}/\code{"mols"} (two names for one estimator; see \sQuote{Details}), which select the modified ordinary least squares moment estimator instead of maximum likelihood: the OLS slopes are kept, the second and third central moments of the OLS residuals are inverted for the scale parameters, and the intercept is shifted up by the implied \eqn{E[u]}. Closed form, no optimizer, deterministic. Implemented for \code{"NHN"}, \code{"NE"} and \code{"NG"}; other models error, because the moment inversion is distribution-specific.

Two further values belong to the wrong-skewness suite (see \link{sfa-wrongskew}) and are implemented for \code{"NHN"} and \code{"NE"} only:

\describe{
  \item{\code{"acols"}}{The alternative corrected OLS of Parmeter and Zhao (2023). As \code{"cols"}, except that the first absolute moment \eqn{E|\varepsilon|} replaces the third central moment. The third moment is the only ingredient that can drive \eqn{\sigma_u} to zero, and it is what fails in small samples. The moment equation has no closed form and is solved numerically. Standard errors need \code{cols_boot}: the analytic errors reported for \code{"cols"} do not apply.}
  \item{\code{"cmle"}}{The moment-constrained maximum likelihood of Zhao and Parmeter (2022): ordinary maximum likelihood subject to two equality constraints the assumed density implies, \eqn{E|\varepsilon| = \sigma\sqrt{2/\pi}} and \eqn{\mathrm{Var}(\varepsilon) = \sigma_v^2 + \frac{\pi-2}{\pi}\sigma_u^2} for \code{"NHN"}, and the analogues of their Proposition 2.1 and Eq. (2.8) for \code{"NE"}. The constraints keep the search away from the Waldman (1982) boundary where OLS is a stationary point of the likelihood. The density is the one \code{estimator = "mle"} maximizes, and the two estimators are asymptotically identical.}
}

See \sQuote{Details}.}
  \item{z_link}{Scale of the variance-determinant linear predictor in the \code{_Z} models. \code{"sd"} (the default, \code{sfm}'s historical convention) sets \eqn{\sigma_u = \exp(z'\delta)}; \code{"var"} sets \eqn{\sigma_u = \sqrt{\exp(z'\delta)}}, as \code{\link{psfm}} and the competing packages do. Both reach the same log-likelihood and the same \code{\link{marginal_effects}}, but \eqn{\delta} differs by a factor of two. Set \code{z_link = "var"} to compare \eqn{\delta} with \code{psfm} or \pkg{sfaR}. See \sQuote{Details}.}
\item{vhet}{One-sided formula for the determinants of the \emph{noise} scale, e.g. \code{vhet = ~ z1 + z2}, giving \eqn{\sigma_v = \exp(z_v'\delta_v)} under the SD link. \code{NULL} (the default) leaves \eqn{\sigma_v} constant. Available for \code{"NHN"}, \code{"NHN_Z"}, \code{"NE"}, \code{"NE_Z"} and \code{"NTN"}.}
\item{uhet}{One-sided formula for the determinants of \eqn{\sigma_u}. The same specification as the \code{| z} segment of \code{formula}; giving both is an error. \code{uhet} lets \code{"NTN"}, which takes no pipe segment, have a heteroskedastic \eqn{\sigma_u}.}
\item{muhet}{One-sided formula for the determinants of the \emph{pre-truncation mean} of \eqn{u}, \eqn{\mu = z_\mu'\delta_\mu}. This is the Battese and Coelli (1995) inefficiency-effects model and requires \code{model_name = "NTN"}, the only family in which \eqn{\mu} exists. \eqn{\mu} carries no link, since it may take either sign.}
  \item{sim_type}{Sequence used to build the simulated-ML draws, for \code{"NLN"} and \code{"NW"} under \code{integration = "simulate"}. This argument and the four after it are ignored under quadrature. \code{"halton"} (default) reproduces previous behaviour. \code{"sobol"} accepts \code{sim_scrambling} 1--3, giving the scrambled sequence of Bhat (2003). \code{"torus"} is the Korobov sequence and accepts \code{sim_prime}. \code{"uniform"} is pseudorandom, a baseline rather than a recommendation. See \sQuote{Details}.}
  \item{antithetics}{Logical. If \code{TRUE}, half the draws are taken from the sequence and the other half created as their mirror images (Hammersley and Morton 1956). Default \code{FALSE}.}
  \item{sim_burn}{Number of leading elements of the sequence discarded. Default \code{NULL}, meaning 1000. The initial elements of different sequences are correlated with one another; discarding them removes the correlation.}
  \item{sim_scrambling}{Integer 0--3, \code{sim_type = "sobol"} only, passed to \code{randtoolbox::sobol}. Non-zero gives a scrambled sequence.}
  \item{sim_prime}{Prime for \code{sim_type = "torus"}. Default \code{NULL} uses the generator's own choice.}
  \item{sim_seed}{Seed for the draws, for reproducibility. The caller's random-number stream is restored afterwards.}
  \item{cols_boot}{Integer. Number of nonparametric bootstrap resamples used for standard errors under \code{estimator = "cols"} or \code{"acols"}. Zero (the default) skips the bootstrap: only the OLS slope standard errors are reported, and the scale parameters and corrected intercept carry \code{NA}. Ignored when \code{estimator = "mle"}.}
  \item{rand.cols}{Integer. Seed for the COLS bootstrap, for reproducibility. The caller's random-number stream is restored afterwards. Default \code{NULL}.}
  \item{integration}{How \code{"NLN"} and \code{"NW"} integrate \eqn{u} out of the composed density. \code{"quadrature"} (the default) is deterministic composite Gauss--Legendre, accurate to about \eqn{10^{-12}} in each observation's log-density. \code{"simulate"} is the simulated maximum likelihood these models used before version 1.3.0, controlled by \code{Nsim} and the \code{sim_*} arguments, and is kept so earlier results can be reproduced. Ignored by every other model. See \sQuote{Details}.}
}
\details{\code{sfm()} fits a cross-sectional frontier \eqn{y = x'\beta + v - u} (production; \eqn{+u} for a cost frontier) by maximum likelihood unless stated otherwise. Choose the distributions of the noise \eqn{v} and the inefficiency \eqn{u \ge 0}:

\tabular{lll}{
 \bold{\code{model_name}} \tab \bold{noise \eqn{v} / inefficiency \eqn{u}} \tab \bold{use it when} \cr
 \code{"NHN"} \tab normal / half-normal \tab the default; the standard frontier \cr
 \code{"NHN_Z"} \tab normal / half-normal, \eqn{\sigma_u = \exp(z'\delta)}{sigma_u = exp(z'delta)} \tab covariates drive the scale of inefficiency \cr
 \code{"NE"} \tab normal / exponential \tab a common alternative to \code{"NHN"} \cr
 \code{"NE_Z"} \tab normal / exponential, \eqn{\sigma_u = \exp(z'\delta)}{sigma_u = exp(z'delta)} \tab as \code{"NHN_Z"}, exponential \eqn{u} \cr
 \code{"NTN"} \tab normal / truncated normal \tab inefficiency has a non-zero mode; with \code{muhet}, Battese-Coelli (1995) \cr
 \code{"NR"} \tab normal / Rayleigh \tab a less skewed alternative to \code{"NHN"} \cr
 \code{"NU"} \tab normal / uniform \tab inefficiency is bounded \cr
 \code{"NG"} \tab normal / gamma \tab a flexible shape for \eqn{u} \cr
 \code{"NNAK"} \tab normal / Nakagami \tab a flexible shape that nests \code{"NR"} \cr
 \code{"NGE"} \tab normal / generalized exponential \tab \eqn{u} may have an interior mode; nests \code{"NE"} \cr
 \code{"TSL"} \tab normal / truncated skew-Laplace \tab as \code{"NGE"}, a different family \cr
 \code{"NLN"} \tab normal / lognormal \tab a lognormal \eqn{u}; computed by quadrature \cr
 \code{"NW"} \tab normal / Weibull \tab a Weibull \eqn{u}; computed by quadrature \cr
 \code{"NGB2"} \tab normal / generalized beta II \tab one fit that nests most of the above \cr
 \code{"tHN"} \tab Student \eqn{t} / half-normal \tab heavy-tailed noise only \cr
 \code{"THT"} \tab Student \eqn{t} / half \eqn{t} \tab heavy tails in both components \cr
 \code{"NB"} \tab normal / binomial \tab wrong skewness read as widespread inefficiency; moments only
}

The subsections below cover the formula and variance determinants, each model family, the other estimators, starting values, boundary fits, standard errors and efficiency prediction. The Monte Carlo measurements behind several of the caveats are collected in \code{vignette("sfa_numerical_notes", package = "sfa")}.

\subsection{Formula and variance determinants}{
The formula is \code{y ~ x | z}. The optional segment after the pipe holds the determinants of \eqn{\sigma_u} for the \code{_Z} models.

\bold{Link scale.} \code{sfm}'s \code{_Z} models place \eqn{z'\delta} on the standard deviation, \eqn{\sigma_u = \exp(z'\delta)}. \code{\link{psfm}}'s \code{TRE_Z}/\code{GTRE_Z}, and the competing packages, place it on the variance, \eqn{\sigma_u = \sqrt{\exp(z'\delta)}}. The two conventions give the same \eqn{\sigma_u}, log-likelihood and marginal effects, but \eqn{\delta} under the SD link is exactly \bold{half} \eqn{\delta} under the variance link. Use \code{z_link} to align them. The default keeps \code{sfm}'s convention so existing results do not move.

\bold{Heteroskedasticity in more than one component.} \code{vhet}, \code{uhet} and \code{muhet} let \eqn{\sigma_v}, \eqn{\sigma_u} and (for \code{"NTN"}) the pre-truncation mean \eqn{\mu} each depend on their own covariates, in any combination. They are named arguments, not further pipe segments, because pipe position already means different things in different families: the second segment is \eqn{\sigma_w} in \code{\link{ttsfm}} but \eqn{\sigma_h} in \code{\link{psfm}}'s \code{GTRE_Z}. The \code{| z} segment keeps its meaning.

Ignored heteroskedasticity in \eqn{v} biases the estimated frontier, not just its precision (Caudill, Ford and Gropper 1995). A \eqn{\sigma_v} that varies with observed characteristics can also absorb apparent wrong skewness without forcing \eqn{\sigma_u} onto zero. \code{muhet} with \code{"NTN"} is the Battese and Coelli (1995) inefficiency-effects model: covariates shift the location of inefficiency, so \eqn{\delta_\mu} has a sign that reads directly. That is usually the model to use for \dQuote{what makes firms inefficient}.

These fits carry every positive quantity on the log scale, so the parameters are unconstrained. They use \code{nlminb} followed by an unbounded BFGS rather than the bounded scaffold the other models use. Coefficients are prefixed by block (\code{Zv.}, \code{Zu.}, \code{Zmu.}), since a covariate may appear in more than one. \code{\link{marginal_effects}} reads the \eqn{\sigma_u} block directly. \code{estimator = "cols"} and the \code{robust} divergences refuse these arguments.

\code{scaling} imposes the scaling property instead, and \code{shapehet} (experimental) lets covariates move the shape of \eqn{u}; see \sQuote{Arguments}.
}

\subsection{Rayleigh, uniform and generalized exponential}{
\bold{\code{"NR"} is not \code{"NHN"}.} The Rayleigh density carries an extra factor of \eqn{u}, and no rescaling removes it. Its standardized skewness is \eqn{-0.631} against the half-normal's \eqn{-0.995}, so data generated under one cannot be fitted by the other. \eqn{\sigma_u} is reported on the second-moment convention \eqn{E[u^2] = \sigma_u^2}, as for \code{"NHN"}, so the two are comparable. The Rayleigh scale is then \eqn{\sigma_u/\sqrt{2}}, with \eqn{E[u] = \sigma_u\sqrt{\pi}/2} and \eqn{Var(u) = (1-\pi/4)\sigma_u^2}. Starting values invert these moments. Test it with \code{data_gen_cs()}'s \code{y_pcs_r} column, not \code{y_pcs}.

\code{"NU"} assumes \eqn{u \sim U(0,\theta)} (Li 1996; Nguyen 2010), giving \eqn{f(e) = \theta^{-1}[\Phi((e+\theta)/\sigma_v) - \Phi(e/\sigma_v)]}. \code{theta} is the upper bound of inefficiency, so predicted inefficiency is bounded by construction.

\code{"NGE"} assumes \eqn{F(u) = (1-e^{-\lambda u})^2}. It nests the exponential model and allows a non-monotonic inefficiency density. \code{sigu} is on the same scale as \code{"NE"}'s (the exponential mean, \eqn{\lambda = 1/\sigma_u}).
}

\subsection{Gamma and Nakagami}{
\code{"NG"} identifies \eqn{E[u] = \mu\sigma_u} far better than \eqn{\mu} or \eqn{\sigma_u} alone (Ritter and Simar 1997). Its start search holds \eqn{E[u]} at a moment estimate and sweeps the shape along that ridge; \code{ng_starts} records the search. The default start can still converge quietly to a much worse optimum. When a fit looks implausible, refit with \code{start_from} set to an \code{"NE"} fit and compare log-likelihoods. \code{"NNAK"}'s shape is weakly identified in the same way.
}

\subsection{Truncated skew-Laplace}{
\code{"TSL"} assumes a truncated skew-Laplace inefficiency (Wang 2012, eq. 3.7; the distribution is due to Aryal and Rao 2005, and the model was published as Nguyen, Gupta and Wang 2014), with density \eqn{f(u) = \frac{1+\lambda}{\sigma_u(2\lambda+1)}\left[2e^{-u/\sigma_u} - e^{-(1+\lambda)u/\sigma_u}\right]} for \eqn{u \ge 0}. This is a signed mixture of two exponentials. It nests the exponential model as \eqn{\lambda \to 0} and allows a non-monotonic density. The composed density is evaluated as a difference in log space, because both terms carry \eqn{\sigma_v^2/(2\sigma_u^2)}, which overflows when \eqn{\sigma_u} is small relative to \eqn{\sigma_v}.

\bold{\code{lambda} is weakly identified, by the family itself.} \eqn{E[U]/\sigma_u} equals 1 at both \eqn{\lambda \to 0} and \eqn{\lambda \to \infty}, and never exceeds \eqn{4-2\sqrt{2} \approx 1.17}; the standardized skewness moves only between about 1.79 and 1.95. The data identify \eqn{\sigma_u \times m(\lambda)} far better than either factor, and the two trade off. With both free, neither \eqn{\sigma_u} nor \eqn{\lambda} attains the root-\eqn{n} rate in simulation; each does when the other is fixed. Report \eqn{E[U]} and \eqn{E[e^{-U}]}, which do converge at close to that rate, rather than the primitives. If \eqn{\lambda} itself matters, use profile-based inference rather than its Wald \eqn{t}-value.
}

\subsection{Lognormal and Weibull}{
\code{"NLN"} (\code{sigu} = sdlog, \code{mu} = meanlog) and \code{"NW"} (\code{sigu} = scale, \code{k} = shape) have no closed-form composed density. By default it is computed by composite Gauss--Legendre quadrature. Each panel is integrated in \eqn{p = F_u(u)} where the inefficiency density is the sharp factor and in \eqn{t = \log u} elsewhere, with breakpoints at the quantiles of \eqn{u}, around the peak of the noise kernel, and around the mode of the integrand. It is accurate to about \eqn{10^{-12}} in log-density wherever that density exceeds \eqn{e^{-50}}, including a Weibull shape below 1, whose density is infinite at zero. The efficiency predictor is the posterior mean over the same nodes.

\code{integration = "simulate"} restores the simulated maximum likelihood of earlier versions. Its log-likelihood error at the default \code{Nsim} is large enough to move \code{sigu} and the shape parameter away from the maximum likelihood estimate, so use it only to reproduce old results.

Under simulation, \eqn{f(e) = \int_0^\infty f_u(u)\,\phi((e+u)/\sigma_v)/\sigma_v\,du} is a product of two densities, and drawing from either alone fails when the other is the narrow factor. Half the draws come from each, combined by the balance heuristic of Veach and Guibas (1995). The draws are generated once and reused at every iteration, so the objective is deterministic. The efficiency predictor reuses the same weights, so \code{u_hat} is the Bayes rule under the density that was maximized.

Each observation gets its own \emph{contiguous block} of the sequence named by \code{sim_type}, not a strided subsequence: a strided van der Corput subsequence is not equidistributed, and blocking is also what induces the negative correlation across observations that makes low-discrepancy draws efficient (Train 2002, ch. 9). Draws are clamped away from 0 and 1, where \code{qnorm()} is infinite.
}

\subsection{Heavy-tailed noise: tHN and THT}{
In \code{"THT"} (Tancredi 2002) one scale mixture is shared by both components, so noise and inefficiency are both Student \eqn{t} with the same degrees of freedom, and the composed error is a closed-form skew-\eqn{t}. In \code{"tHN"} the noise is Student \eqn{t} and the inefficiency is half-normal, drawn \emph{independently}. There is no closed form; the density \eqn{f(e) = \int_0^\infty f_v(e+u) f_u(u) du} is evaluated by Gauss--Legendre quadrature. \code{"tHN"} is the heavy-tailed-\emph{noise} model, and so the natural parametric comparison for the robust estimators; \code{"THT"} is not, because its inefficiency is heavy-tailed too. \code{"tHN"} reports \code{(sigv, sigu, nu)}, not \code{"THT"}'s inverted \code{(sigu, sigv, a)}.

Caveats for \code{"tHN"}:
\itemize{
  \item The degrees of freedom are weakly identified: the profile likelihood is nearly flat over a wide range of \eqn{\nu}. Do not take a single \eqn{\nu} at face value; profile over a grid of fixed \eqn{\nu} and report the path. \code{sfm()} fits from several separated starts, keeps the best, records them in \code{thn_starts}, and warns when they reach different optima.
  \item On real data \eqn{\sigma_u} can collapse onto zero, the heavy noise tail absorbing the one-sided component. That is a property of the model, so \code{sfm()} warns and sets \code{thn_sigma_u_at_bound}.
  \item Fits are slower than the closed-form models: every evaluation is a quadrature whose node count grows with \eqn{\sigma_u/\sigma_v}.
}
}

\subsection{Generalized beta of the second kind (NGB2)}{
Makiela and Mazur (2022), section 3, with density
\deqn{f(u) = \frac{\psi_u}{\sigma_u B(\tau/\psi_u, \nu_u/\psi_u)\,\nu_u^{\tau/\psi_u}}
      \left(\frac{u}{\sigma_u}\right)^{\tau-1}
      \left[1 + \frac{1}{\nu_u}\left(\frac{u}{\sigma_u}\right)^{\psi_u}\right]^{-(\tau+\nu_u)/\psi_u},}
reported as \code{(sigv, sigu, nu, psi, tau)}. \eqn{\sigma_u} is a scale,
\eqn{\nu_u} governs tail thickness (and so which moments exist), and
\eqn{\psi_u} and \eqn{\tau} are shape. Its value is what it nests, so one fit
says which of the simpler families the data want:

\tabular{ll}{
 \bold{restriction} \tab \bold{distribution} \cr
 \eqn{\nu_u \to \infty} \tab generalized gamma \cr
 \eqn{\nu_u \to \infty}, \eqn{\tau = 1}, \eqn{\psi_u = 2} \tab half-normal (\code{"NHN"}) \cr
 \eqn{\nu_u \to \infty}, \eqn{\tau = \psi_u = 1} \tab exponential (\code{"NE"}) \cr
 \eqn{\nu_u \to \infty}, \eqn{\psi_u = 1} \tab gamma (\code{"NG"}) \cr
 \eqn{\nu_u \to \infty}, \eqn{\tau = \psi_u} \tab Weibull (\code{"NW"}), at scale \eqn{\sigma_u\psi_u^{1/\psi_u}} \cr
 \eqn{\tau = 1}, \eqn{\psi_u = 2} \tab half Student \eqn{t} with \eqn{\nu_u} degrees of freedom \cr
 \eqn{\tau = \psi_u = 1} \tab Lomax (Pareto II), shape \eqn{\nu_u} and scale \eqn{\sigma_u\nu_u}
}

The half-\eqn{t} and Lomax cases are exact; the others are limits. The Lomax
row is the exponential row with \eqn{\nu_u} left free, as the half-\eqn{t} row
is the half-normal row with \eqn{\nu_u} free, so \code{"NGB2"} already covers a
heavy-tailed Lomax model.

Caveats:
\itemize{
  \item \bold{The Lomax region can lack a mean.} A Lomax with
  \eqn{\nu_u \le 1} has no finite mean, and with \eqn{\nu_u \le 2} no finite
  variance. A fit near \eqn{\tau = \psi_u = 1} with small \eqn{\nu_u} is
  well defined, but \code{exp_u_hat} and the JLMS predictor then estimate a
  population mean that does not exist.
  \item \bold{A large \eqn{\nu_u} is not an estimate.} Once the tail is not
  heavy the likelihood is flat in \eqn{\nu_u}, with no interior maximum.
  \code{sfm()} warns above 200. Read that as \dQuote{generalized gamma, no
  heavy tail detected}, treat the standard error on \eqn{\nu_u} as
  meaningless, and read \eqn{\psi_u} and \eqn{\tau}. Starts are at the
  half-normal point (\eqn{\tau = 1}, \eqn{\psi_u = 2}, \eqn{\nu_u = 30}), so
  a fit that moves away from it has moved for a reason.
  \item \bold{It is the slowest \code{sfm()} model}: a fit takes minutes, not seconds.
}

The composed density is evaluated by deterministic quadrature with two node
sets, one at the quantiles of \eqn{u} and one at the quantiles of the noise
truncated to \eqn{u > 0}, combined by the balance heuristic. Either alone
fails: the first misses a narrow noise, the second a concentrated \eqn{u}.
It is deliberately \emph{not} simulated: with four shape parameters, a
simulated likelihood lets the optimizer move toward parameters that maximize
its own simulation error.
}

\subsection{Binomial inefficiency (NB)}{
Carree (2002). Every other one-sided distribution here is positively skewed,
so a production frontier implies a negatively skewed composed error, and a
positive residual skew leads to \eqn{\sigma_u = 0} (Waldman 1982), read as
\dQuote{no inefficiency}. A binomial \eqn{b(n, p)} is negatively skewed for
\eqn{p > 1/2}, which says the opposite: most firms carry substantial
inefficiency and only a few sit near the frontier.

The estimator inverts the second, third and fourth central moments of the OLS
residuals (his Eq. 7),
\deqn{m_2 = \sigma_v^2 + np(1-p), \quad m_3 = -np(1-p)(1-2p), \quad
      m_4 - 3m_2^2 = np(1-p)(1-6p+6p^2),}
whose ratio \eqn{x = (m_4 - 3m_2^2)/m_3 = (6p^2-6p+1)/(2p-1)} is a quadratic in
\eqn{p} with roots \eqn{p = \tfrac12 + \tfrac{x}{6} \pm \tfrac16\sqrt{x^2+3}}.
The root is chosen by \eqn{x} outside \eqn{[-1,1]} and by the sign of
\eqn{m_3} inside it. \code{n_bin} and \code{p_bin} are reported with
\code{sigu} \eqn{= \sqrt{np(1-p)}}, comparable with the other models. The
efficiency predictor is the exact discrete posterior
\eqn{P(u = j \mid \varepsilon) \propto \binom{n}{j}p^j
(1-p)^{n-j}\phi(\varepsilon + j; 0, \sigma_v)}, with \eqn{n} rounded to an
integer; \code{n_bin} keeps the unrounded estimate.

Caveats:
\itemize{
  \item \bold{A positive \code{m3} is admissible}; that is the point of the
  model. \code{wrong_skew} is set only when the inversion fails.
  \item \bold{It uses the fourth moment}, so it is more fragile than the
  half-normal and exponential inversions. In Carree's infeasible region
  \eqn{|m_4 - 3m_2^2| > |m_3| > 0} the implied \eqn{n} is negative;
  \code{sfm()} warns and reports no estimate. Read that as evidence against a
  binomial inefficiency. It also needs \eqn{\sigma_v} not small relative to
  \eqn{\sqrt{np(1-p)}}, since \eqn{\sigma_v^2} is then a small difference of
  large moments.
  \item \bold{The support is the integers}: one unit of \eqn{u} is one unit of
  the dependent variable, so with a log dependent variable each category is a
  factor of \eqn{e}. No scale factor is identified, so read the implied
  efficiency levels with that in mind.
}
}

\subsection{Estimators other than maximum likelihood}{
\bold{Modified OLS (\code{estimator = "cols"} or \code{"mols"}).} The moment estimator of Olson, Schmidt and Waldman (1980), usually called MOLS. It is \emph{not} Winsten's corrected OLS, which shifts the intercept by the largest residual and estimates no variance parameters. The package has always called it \code{"cols"}, which still works; \code{"mols"} is an exact synonym. OLS is consistent for the slopes whatever the one-sided distribution; only the intercept is biased, by \eqn{E[u]}. MOLS keeps the OLS slopes, inverts the central moments of the residuals for the scale parameters, and shifts the intercept by the implied \eqn{E[u]}. With \eqn{m_2}, \eqn{m_3} the second and third central moments of \eqn{e = v - u}: for \code{"NHN"}, \eqn{m_3 = \sigma_u^3\sqrt{2/\pi}(1-4/\pi)} and \eqn{m_2 = \sigma_v^2 + \sigma_u^2(1-2/\pi)} with \eqn{E[u] = \sigma_u\sqrt{2/\pi}}; for \code{"NE"}, \eqn{m_3 = -2\sigma_u^3} and \eqn{m_2 = \sigma_v^2+\sigma_u^2} with \eqn{E[u] = \sigma_u}; for \code{"NG"}, from the third and fourth cumulants, \eqn{\sigma_u = -k_4/(3k_3)} and \eqn{\mu = -k_3/(2\sigma_u^3)} with \eqn{E[u] = \mu\sigma_u}.

No optimizer runs and the answer is deterministic, so MOLS is a useful check on a maximum-likelihood fit that may sit at a local optimum. It is less efficient than ML when the distribution is right, and the \code{"NG"} inversion, which uses a fourth cumulant, is noisier than the other two at moderate sample sizes.

\bold{Wrong skew under MOLS.} A production frontier implies \eqn{m_3 < 0}. When \eqn{m_3 \ge 0} the moment equations have no admissible solution: \code{sfm()} warns, reports \code{sigu} as 0 with all residual variance in \code{sigv}, and returns no efficiency predictions. Read that as \emph{no evidence of inefficiency in these data}, not an estimate of zero. This is the Type I failure of Olson, Schmidt and Waldman, common in small samples; \code{$wrong_skew} and \code{$residual_moments} record it.

\bold{The wrong-skewness suite (\code{"acols"} and \code{"cmle"}).} Waldman (1982) showed that when OLS residuals are skewed the wrong way, the OLS fit is a local maximum of the likelihood at \eqn{\hat\sigma_u = 0}. Wrong skew arises routinely when \eqn{\sigma_u} is small relative to \eqn{\sigma_v}, within the assumed model (Simar and Wilson 2009). \code{"acols"} removes the third moment from the moment estimator, the only ingredient that can produce the boundary. \code{"cmle"} keeps the likelihood and adds two constraints the assumed density implies, which the boundary does not satisfy. Neither changes the model. See \link{sfa-wrongskew} for the other options.
}

\subsection{Robust estimation}{
For \code{model_name = "NHN"}, \code{robust} replaces the log-likelihood with
a robust divergence objective that downweights outlying observations:
\itemize{
  \item \code{"mlqe"}: minimum L\eqn{q} estimation (Ferrari and Yang 2010),
  power parameter \code{c_mlqe} (\eqn{q = 1-c}). Objective
  \eqn{\sum (f(e_i)^c - 1)/c}.
  \item \code{"psi"}: Psi-likelihood, power parameter \code{eta}. Objective
  \eqn{\sum f(e_i)^c/c - n \, I_{1+c}/(1+c)}, where
  \eqn{I_{1+c} = \int f(e;\theta)^{1+c}\,de} is evaluated numerically at
  every iteration.
  \item \code{"mdpd"}: minimum density power divergence (Basu et al. 1998),
  power parameter \code{alpha}. Objective
  \eqn{((1+c)/c) \sum f(e_i)^c - n \, I_{1+c}}.
}
All three reduce to \code{"mle"} as the power parameter goes to 0.
\strong{Psi and MDPD give identical point estimates for the same power
parameter} (\eqn{c = \eta = \alpha}): \eqn{\text{MDPD}(c) = (1+c)\cdot
\text{Psi}(c)}, a positive rescaling. Both are offered to match the source
paper's terminology. MLqE is a distinct estimator.

These are M-estimators, so the Hessian-inverse standard error is not valid
for them. \code{sfm()} reports the sandwich form \eqn{A^{-1} B A^{-1}}
instead: \eqn{A} is the Hessian of the summed objective and
\eqn{B = \sum_i g_i g_i'} is built from per-observation gradients
(\code{\link[numDeriv]{jacobian}}). This requires \code{optHessian = TRUE}.
}

\subsection{Starting values and optimization}{
Most models run a staged optimizer: \code{nlminb} or bobyqa (see \code{use.nlminb}, \code{use.bobyqa}), then optionally \code{psoptim} (\code{PSopt}), then \code{optim} with \code{"L-BFGS-B"}, each stage seeding the next. The iteration defaults are high because \code{sfm()}'s likelihoods are cheap to evaluate. \code{\link{psfm}}'s are lower because several of its models use simulated likelihoods.

Supply \code{start_val} to fix the start, or \code{start_from} to seed from another fit by parameter name. Cascading from a simpler model is the main defence against a fit that stops at a worse optimum without warning, most often for \code{"NG"}.

After the fit, \code{sfm()} checks that the returned point is a maximum (\code{stationarity} under \sQuote{Value}). \code{optim()} can stop with code 0 on a flat ridge, so when the Newton decrement is too large \code{sfm()} polishes with \code{nlminb()} and keeps the result only if the log-likelihood rises.
}

\subsection{Boundary fits and wrong skew}{
A scale parameter on its zero boundary is often the maximum likelihood estimate for that sample, not a numerical failure. \code{sfm()} reports it rather than bounding the parameter away from zero:
\itemize{
  \item \code{sigma_u_at_bound} with \code{wrong_skew}: the residuals are skewed the wrong way, and \eqn{\sigma_u = 0} is the estimate (Waldman 1982). Read it as no evidence of inefficiency. With correctly skewed residuals, a collapse is suspect and \code{sfm()} warns.
  \item \code{thn_sigma_u_at_bound} (\code{"tHN"}): the noise tail absorbed the inefficiency.
  \item \code{ntn_sigma_v_at_bound} (\code{"NTN"}) and \code{sigma_v_at_bound} (\code{"NHN_Z"}, \code{"NE_Z"}, \code{"NR"}): \eqn{\sigma_v} collapsed and the fit is a deterministic frontier. Standard errors are \code{NA}.
}
Each is described under \sQuote{Value}. In every case treat the efficiency scores as uninformative. See \link{sfa-wrongskew} and \code{\link{skewness_test}} for tests and remedies.
}

\subsection{Standard errors}{
Standard errors come from the numerical Hessian and need \code{optHessian = TRUE}; otherwise they are \code{NA}. Under \code{robust != "mle"} they are sandwich-form (see \sQuote{Robust estimation}). Under sampling \code{weights} use \code{sandwich::sandwich(fit)}. Under \code{estimator = "cols"} the scale parameters need \code{cols_boot > 0}; \code{"acols"} always needs it. At a boundary fit no standard error exists and \code{NA} is reported.
}

\subsection{Efficiency prediction}{
\code{exp_u_hat} is \eqn{E[\exp(-u)|e]} (Battese and Coelli 1988) and \code{u_hat} is the Jondrow et al. (1982) predictor \eqn{E[u|e]}, each under the posterior the model's own \eqn{u} implies. \sQuote{Value} lists which models return each. \code{"THT"} also returns \code{sd_exp_u_hat}, and \code{"NHN"} returns \code{med_u_hat}. See \code{\link{efficiency}} for a common interface.
}
}
\author{David H. Bernstein and Alexander Stead}
\note{Standard errors require optHessian set to TRUE. When
\code{robust != "mle"} (\code{model_name = "NHN"} only), a sandwich-form SE
is reported instead of the ordinary Hessian-inverse one; see Details.}

\references{
Parmeter, C. F. and Zhao, S. (2023). An alternative corrected ordinary least squares estimator for the stochastic frontier model. \emph{Empirical Economics} 64(6), 2831-2857. \doi{10.1007/s00181-023-02401-1}

Zhao, S. and Parmeter, C. F. (2022). The \dQuote{wrong skewness} problem: Moment constrained maximum likelihood estimation of the stochastic frontier model. \emph{Economics Letters} 221, 110901. \doi{10.1016/j.econlet.2022.110901}

Waldman, D. M. (1982). A stationary point for the stochastic frontier likelihood. \emph{Journal of Econometrics} 18(2), 275-279. \doi{10.1016/0304-4076(82)90041-0}

Simar, L. and Wilson, P. W. (2009). Inferences from cross-sectional, stochastic frontier models. \emph{Econometric Reviews} 29(1), 62-98. \doi{10.1080/07474930903324523}
}

\seealso{\link{sfa-wrongskew} for the full wrong-skewness suite, \code{\link{psfm}} for panel models, \code{\link{zsfm}} for zero-inefficiency models, \code{\link{ttsfm}} for two-tier models, \code{\link{data_gen_cs}} for simulating data with known true parameters, \code{\link{sfareg-methods}} for extracting results, and \code{vignette("sfa_numerical_notes", package = "sfa")} for the measurements behind the caveats above.}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{exp_u_hat}{Predicted technical efficiency, E[exp(-u)|e]. Available for models: NHN, NHN_Z, NR, NG, NNAK, NE, NTN, NU, NGE, NLN, NW, NGB2, THT and tHN, and for \code{estimator = "cols"} with NHN, NE, NG and NB. Each uses the posterior its own \eqn{u} implies.}
\item{u_hat}{Predicted inefficiency, E[u|e] (Jondrow et al. 1982). Available for NE, NTN, NU, NGE, NLN, NW, NGB2, THT and tHN.}
\item{n_bin, p_bin}{\code{model_name = "NB"} only, in \code{out} and \code{coefficients}: the moment estimates of the binomial's size and success probability. \code{n_bin} is a real number, not an integer.}
\item{wrong_skew, sigma_u_at_bound, residual_m3}{Wrong-skewness diagnostics. \code{residual_m3} is the third central moment of the OLS residuals. \code{wrong_skew} is \code{TRUE} when it is non-negative, which for a production frontier means the moment equations have no interior solution. \code{sigma_u_at_bound} is \code{TRUE} when the fitted one-sided scale has collapsed to zero. When both are \code{TRUE} the boundary \emph{is} the maximum likelihood estimate (the Type I failure of Olson, Schmidt and Waldman 1980; Waldman 1982), so \code{sfm} reports it rather than bounding \eqn{\sigma_u} away from zero. Read it as no evidence of inefficiency, and treat the efficiency scores as uninformative. A collapse with correctly skewed residuals is the opposite case: the boundary is unlikely to be the maximum, and \code{sfm} warns that the fit is suspect.}
\item{thn_sigma_u_at_bound}{Logical, \code{tHN} only. \code{TRUE} when the fitted \eqn{\sigma_u} collapsed onto zero, so the model reports essentially no inefficiency. Accompanied by a warning. Treat the efficiency scores as uninformative when this is set.}
\item{ntn_sigma_v_at_bound}{Logical, \code{NTN} only. \code{TRUE} when the fitted \eqn{\sigma_v} collapsed onto zero (\code{lambda} ran away), so the fit attributes all variation to inefficiency. On such a sample the likelihood rises monotonically toward a finite limit as \code{lambda} grows, so the boundary is the maximum likelihood estimate. The Hessian is singular there and every standard error is \code{NA}. Accompanied by a warning. Read it as \dQuote{no noise component identified}, not as an estimate of \code{lambda}. Because the optimizer can stop at an interior \code{lambda} below the supremum, the fit climbs a ladder of fixed \code{lambda} values from its own solution and moves to the corner when the likelihood there is higher by more than 1e-4.}
\item{sigma_v_at_bound}{Logical, \code{NHN_Z}, \code{NE_Z} and \code{NR} only. \code{TRUE} when the fitted \eqn{\sigma_v} collapsed onto its lower bound relative to the median \eqn{\sigma_u}, so the fit is a deterministic frontier. The optimizer found no interior point with a higher likelihood, so the boundary is the maximum likelihood estimate on that sample. Every standard error is \code{NA}. For \code{NR} the fit also checks the corner when the optimizer stops inside: it maximizes with \eqn{\sigma_v} fixed near zero and moves there only if the likelihood is higher by more than 1e-4. Accompanied by a warning.}

\item{se_at_bound}{Character, present for a maximum likelihood fit only when a parameter lies within \code{optim()}'s Hessian differencing step (2e-3) of its lower bound, as \code{"NG"} and \code{"NNAK"} scale and shape parameters often do. Those parameters have no standard error, because the Hessian would difference across the bound; the other standard errors come from the Hessian over the remaining parameters and are conditional on the bound ones.}
\item{se_unresolved}{Present, for a maximum likelihood fit with \code{optHessian = TRUE}, only when at least one reported standard error could not be computed and no boundary flag (\code{sigma_u_at_bound}, \code{sigma_v_at_bound}, \code{ntn_sigma_v_at_bound}, \code{thn_sigma_u_at_bound}) already explains it. A list naming those parameters (\code{parameters}, \code{n_missing} of \code{n_par}), the Hessian's \code{min_eigenvalue} and \code{max_eigenvalue}, whether it is \code{indefinite}, and the \code{message} the fit warned with. An indefinite Hessian means the estimated variance is negative, so no standard error exists at that point. It is either a point that is not a maximum or a boundary solution, and this component does not distinguish the two. Printed by \code{print()} and \code{summary()}, because a warning does not survive \code{saveRDS()}.}
\item{stationarity}{List, present when standard errors were to be reported from a maximum likelihood fit. \code{decrement} is the Newton decrement \eqn{\frac{1}{2} g^\top H^{-1} g} at the returned point (the log-likelihood gain one Newton step predicts), from a Richardson gradient; \code{grad_max} is the largest absolute gradient element. A coordinate closer to a bound than the differencing stencil is differenced forward instead, left out of the decrement, and judged by the sign of its one-sided slope. \code{at_bound} marks those coordinates, \code{n_at_bound} counts them, and \code{bound_ok} is \code{TRUE} when each has the non-negative slope a constrained minimum of the negative log-likelihood requires. A \emph{negative} decrement means the Hessian over the remaining coordinates is not positive definite, so the point is a saddle. When the decrement exceeds \code{1e-4} in absolute value, \code{sfm} polishes with \code{nlminb()}, positive parameters on the log scale, and keeps the result only if the log-likelihood rises; \code{polished}, \code{gain} and \code{decrement_after} record what happened. A polish that reaches a boundary triggers that boundary's flag and warning. If the point is still not a maximum and no boundary explains it, \code{sfm} warns.}

\item{thn_starts}{List, \code{tHN} only. Multi-start diagnostic: how many starting values were tried and converged, the log-likelihood reached from each, and \code{n_distinct}, the number of separate optima found. \code{n_distinct > 1} means the likelihood is multimodal for these data and the returned fit is only the best of those found.}
\item{estimator}{\code{"cols"}, \code{"acols"} or \code{"cmle"}, naming the estimator that ran; absent for ordinary maximum-likelihood fits.}
\item{wrong_skew}{Logical, \code{estimator = "cols"} and \code{"acols"} only. For \code{"cols"}, \code{TRUE} when the OLS residuals were skewed the wrong way and the moment equations had no admissible solution. For \code{"acols"} the third moment is never consulted, so it is \code{TRUE} when the absolute-moment equation has no root on the admissible interval.}

\item{boundary}{Logical, \code{estimator = "cmle"} only. \code{TRUE} when the two moment constraints could not both be met in the interior and one of \eqn{\sigma_u}, \eqn{\sigma_v} sits on its boundary.}
\item{residual_moments}{The second, third and fourth central moments of the OLS residuals, \code{estimator = "cols"} only.}
\item{cols_boot_draws}{The matrix of bootstrap replicates (\code{cols_boot} rows, one column per parameter) when \code{cols_boot > 0}; \code{NULL} otherwise.}
\item{ng_starts}{Returned by \code{NG} only. A record of the starting-value search: how many candidates were built (\code{n_tried}), how many were polished before choosing (\code{n_polished}), the log-likelihood at each candidate (\code{loglik_at_start}), and the best value reached (\code{best}). \code{NULL} when \code{start_val} was supplied. The candidates hold \eqn{E[u] = \mu\sigma_u} at a moment estimate and sweep the shape along that ridge (Ritter and Simar 1997). Widely differing log-likelihoods among the candidates mean the surface is more awkward than usual for that data set.}
\item{sd_exp_u_hat}{Standard deviation of exp(-u) given the residual, following Tancredi (2002, section 2.2). Returned by \code{THT} only. For a large \emph{positive} residual the half-normal model concentrates exp(-u) at 1 with near-zero uncertainty, while the skew-t model treats the point as an outlier, widens this spread and does not push the efficiency estimate to 1.}
\item{med_u_hat}{Predicted technical efficiency (median values). Available only for the NHN model.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}). Sandwich-form (not the ordinary Hessian-inverse form) when \code{robust != "mle"}; see Details.}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{data}{The \code{data} argument as supplied, so that \code{fitted()}, \code{residuals()}, \code{predict()} and \code{efficiency(logDepVar = FALSE)} need not find it again, which fails once it is out of scope (for example, a fit made inside a function from a formula built outside it). \code{NULL} when no \code{data} was given.}
\item{inefdec}{The orientation this fit used, as a scalar logical: \code{TRUE} for a production frontier, \code{FALSE} for a cost one. Recorded at fit time so that post-estimation functions (\code{\link{efficiency}} with \code{logDepVar = FALSE}, the robust accessors, the composed CDF) need not re-evaluate the \code{inefdec} argument of the recorded call, which fails when it was passed as a variable no longer in scope.}
\item{het}{Only present when \code{vhet}, \code{uhet} or \code{muhet} was given: the family fitted, the link used, the sizes of the four parameter blocks, the formulas as supplied, and (for \code{"NTN"}) the fitted per-observation \eqn{\mu}.}
\item{sigma_u, sigma_v}{Only present for a heteroskedastic fit: the fitted per-observation scales.}
\item{v_spec, mu_spec}{Only present for a heteroskedastic fit: the noise and pre-truncation-mean blocks, in the same shape as \code{z_spec}.}
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

## Heteroskedasticity in the noise as well as in the inefficiency: `| z`
## still drives sigma_u, `vhet` drives sigma_v. Coefficients come back
## prefixed by block (Zu., Zv.).
cs.het      <-  sfm(formula = y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z",
                    data = cs_data_trial, vhet = ~ z)

## Battese and Coelli (1995): covariates shift the LOCATION of inefficiency
## rather than its dispersion. Requires the truncated-normal family.
cs.bc95     <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NTN",
                    data = cs_data_trial, muhet = ~ z)
}
}
