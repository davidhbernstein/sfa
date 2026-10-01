\name{ttsfm}
\alias{ttsfm}
\title{Two-Tier Stochastic Frontier Model}
\description{
Fits the two-tier stochastic frontier model (2TSF): the composed error is a symmetric noise term plus two one-sided components pulling in opposite directions, \eqn{\varepsilon = v - u + w} (for example buyer and seller informational rents).
\itemize{
\item \code{"TTNE"}: \eqn{u} and \eqn{w} exponential.
\item \code{"TTHN"}: \eqn{u} and \eqn{w} half-normal.
\item \code{"TTNLS"}: nonlinear least squares, with no distributional assumption. It estimates the frontier slopes only; see \sQuote{Details}.
}
}
\usage{ttsfm(formula, model_name = c("TTNE", "TTHN", "TTNLS"),
data, z_link = c("sd", "var"), maxit.bobyqa = 80000,
maxit.psoptim = 1000, maxit.optim = 1000,
REPORT = 1, trace = 0, pgtol = 0, start_val = FALSE, PSopt = FALSE,
optHessian = TRUE, inefdec = TRUE, upper = NA,
Method = "L-BFGS-B", logit = TRUE, verbose = FALSE,
keep_objective = FALSE, rand.psoptim = NULL)}
\arguments{
\item{formula}{A symbolic description of the frontier, \code{y ~ x1 + x2}. Use \code{y ~ x | zu | zw} to allow determinants \code{zu} of the first one-sided component and \code{zw} of the second.}
\item{model_name}{\code{"TTNE"}, \code{"TTHN"} or \code{"TTNLS"}; see \sQuote{Description}. Matching ignores case, so \code{"ttne"} and \code{"TTNE"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{z_link}{Scale on which the variance-determinant linear predictors are placed. \code{"sd"} (the default, and \code{ttsfm}'s historical convention) sets \eqn{\sigma_u = \exp(z_u'\delta)} and likewise for \eqn{\sigma_w}; \code{"var"} sets \eqn{\sigma_u = \sqrt{\exp(z_u'\delta)}}, matching \code{\link{psfm}} and the competing packages. The two are reparameterizations of the same model -- identical maximised log-likelihood and identical frontier coefficients -- but \eqn{\delta} under the SD link is exactly \bold{half} \eqn{\delta} under the variance link, so a \eqn{\delta} from \code{ttsfm} is not comparable with one from \code{psfm} or \pkg{sfaR} unless the same convention is used. See \code{\link{sfm}}, which gained the same argument.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine.}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine.}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine.}
\item{REPORT, trace, pgtol, upper}{Accepted but not used by \code{ttsfm()}: the optimizer's own tracing is controlled by \code{verbose}, and its bounds are set internally.}
\item{start_val}{\code{FALSE} (the default) to build starting values from the data, or a numeric vector giving them directly, in the optimizer's parameter order.}
\item{PSopt}{Logical. Run the particle-swarm stage (\code{pso::psoptim}) between \code{bobyqa} and \code{optim}. \code{FALSE} by default. The swarm is the only stochastic step in \code{ttsfm()}: with the default, a fit is fully deterministic given the data and does not depend on the session's random-number state. See \code{rand.psoptim}.}
\item{optHessian}{Logical. Finish with \code{optim()} and a numerically differentiated Hessian? Standard errors and t-values require \code{TRUE}, the default.}
\item{inefdec}{Has no effect in \code{ttsfm()}: the composed error is \eqn{v - u + w} whichever value is given, and \code{TRUE} and \code{FALSE} return identical fits. Read \code{sigma_u} as the component that lowers \eqn{y} and \code{sigma_w} as the one that raises it.}
\item{Method}{Method passed to \code{optim} for the final stage. Must be
\code{"L-BFGS-B"}: it is the only \code{optim} method that applies the lower
and upper bounds this package uses to keep scale parameters positive, and any
other value is refused with an error rather than quietly optimized without
them.}
\item{logit}{Accepted but not used by \code{ttsfm()}.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{keep_objective}{Logical. Retain the likelihood closure on the fitted object as \code{$objective}, so the score-based tools can difference it: \code{vcov(type = "bhhh")} (OPG/BHHH standard errors), \code{\link{influence_sfa}}, the \pkg{sandwich} methods, \code{\link{TIC}} and \code{vuong}. New in 1.2.1, for \code{"TTNE"} and \code{"TTHN"}. It has no effect for \code{"TTNLS"}, which is nonlinear least squares rather than maximum likelihood and so has no score; that case warns. \code{FALSE} by default, because a closure carries its enclosing environment.}
\item{rand.psoptim}{Integer seed for the particle-swarm stage, used only when \code{PSopt = TRUE}. The seed is applied to that stage alone and the caller's random-number state is restored afterwards. With \code{PSopt = TRUE} and the default \code{NULL} the swarm draws from the session's random-number state, so set a seed whenever a result must be reproducible across sessions or package versions.}
}
\details{
\subsection{The model}{
The two-tier stochastic frontier model (Polachek and Yoon, 1987) writes the composed error as \eqn{\varepsilon = v - u + w}, with \eqn{v} symmetric normal noise and \eqn{u}, \eqn{w} independent one-sided components: \eqn{u} lowers \eqn{y} below the frontier and \eqn{w} raises it above. \code{"TTNE"} assumes \eqn{u} and \eqn{w} are exponential; \code{"TTHN"} assumes they are half-normal.

Determinants of \eqn{u} and \eqn{w} enter through the second and third parts of a two-pipe formula, \code{y ~ x | zu | zw}, which parameterize \eqn{\sigma_u} and \eqn{\sigma_w} through a log link (see \code{z_link}).
}

\subsection{\code{"TTNLS"}}{
Nonlinear least squares makes no distributional assumption, and the price is that it identifies only the frontier slopes and a composite intercept. \eqn{E[y \mid x] = x'\beta - E[u] + E[w]}, so the sum of squares cannot separate \eqn{\sigma_u}, \eqn{\sigma_w} and the intercept. \code{ttsfm()} therefore:
\itemize{
\item reports an intercept that absorbs \eqn{E[w] - E[u]}. Its standard error is not meaningful: the Hessian it comes from still contains the unidentified \eqn{\sigma_u} and \eqn{\sigma_w} directions (12.9 against 0.013 for the slope, on 3,000 simulated observations);
\item returns \code{NA} for \code{sigma_v}, \code{sigma_u} and \code{sigma_w}, and for any determinant coefficients, with a warning;
\item returns \code{metrics} whose entries are all \code{NA}, since they are functions of those scales.
}
Use \code{"TTNE"} or \code{"TTHN"} when the individual scales are needed.
}

\subsection{Boundary fits and \code{metrics}}{
Like the other composed-error models in this package, \code{"TTNE"} and \code{"TTHN"} can have a genuine (non-error) boundary mode where a variance component is estimated near its numerical floor in some samples; spot-check fits before serious use. \code{"TTHN"}'s information-deficiency metrics (the \code{metrics} component) involve a large number of bivariate-normal-CDF evaluations and are more experimental than the parameter estimates themselves; they are \code{NULL} if that calculation fails.
}

\subsection{Iteration defaults}{
Like \code{sfm()} and \code{zsfm()}, \code{ttsfm()}'s likelihoods are closed-form, so its optimization defaults are much higher than \code{psfm()}'s, which uses simulated maximum likelihood for several models. See \code{\link{psfm}}'s Details.
}
}
\references{Polachek, S. and Yoon, B.J. (1987). A Two-Tiered Earnings Frontier Estimation of Employer and Employee Information in the Labor Market. \emph{The Review of Economics and Statistics}.}
\author{David H. Bernstein and Christopher F. Parmeter}
\note{Standard errors require \code{optHessian = TRUE}.

\bold{\code{"TTHN"} will drive \eqn{\sigma_v} to zero when \eqn{\sigma_u} and \eqn{\sigma_w} are close, and that is the maximum likelihood estimate rather than a failure of the optimizer.} With \eqn{\sigma_u = \sigma_w} the difference \eqn{w - u} is already symmetric and heavy-tailed, so a small normal \eqn{v} underneath it is very nearly invisible: on a 2,000-observation draw at \eqn{\sigma_u = \sigma_w = 1} the \eqn{\sigma_v} profile is flat to 0.000 log units from \eqn{10^{-6}} upward and decreasing above 0.01, with no interior maximum, and the collapsed fit beats the true parameter values by 1.29 log units. The two-tier structure is identified; the noise underneath two equal and opposite tiers is not. Asymmetric \eqn{(\sigma_u, \sigma_w)} does not show this.}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values, one row per parameter.

\bold{Scale parameters are reported on their natural scale and named, as of 1.2.1.} The rows are \code{sigma_v}, \code{sigma_u} and \code{sigma_w}. Before 1.2.1 this matrix held raw optimizer values -- the row labelled \code{sigv} was \eqn{\log\sigma_v}, and in the homoskedastic case the two one-sided scales were \emph{both} labelled \code{(Intercept)}, so \eqn{u} and \eqn{w} could not be told apart. Code that exponentiated these rows must stop doing so.

Where \eqn{\sigma_u} or \eqn{\sigma_w} carries determinants there is no single scale to report, so those rows stay as \emph{coefficients} on the link and take the prefixes \code{Zu.} and \code{Zw.} -- which is also what makes the two blocks distinguishable when they share a covariate name. \code{sigma_v} is homoskedastic in every case and is always reported directly.

\bold{A t-value on a scale parameter is not a conventional test.} The null it corresponds to sits on the boundary of the parameter space, where the usual normal reference distribution does not apply. The column is reported for consistency with the rest of the package; read it as a rough magnitude, not as a p-value.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{\code{"TTNE"}, \code{"TTHN"} or \code{"TTNLS"}.}
\item{formula}{The formula used in the model specification.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{metrics}{A list of information-deficiency measures (M1, M2, M5, M6, M7, M10; see Polachek and Yoon) evaluated at the final parameter estimates. \code{NULL} if the underlying calculation failed, and all \code{NA} for \code{"TTNLS"}; see \sQuote{Details}.}
\item{call}{The matched call.}
\item{nobs}{The number of observations used in the fit.}
\item{objective}{The likelihood closure, when \code{keep_objective = TRUE} (\code{"TTNE"} and \code{"TTHN"} only).}
}
\examples{
\donttest{
library(sfa)

data <- data_gen_cs(N = 200, rand = 1, sig_u = 0.3, sig_v = 0.2,
                     cons = 2, beta1 = 0.5, beta2 = 0.3, a = 5, mu = 0)

eq <- y_pcs ~ x1 + x2

ttsfm(formula    = eq,
      model_name = "TTNE",
      data       = data)
}
}
