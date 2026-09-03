## Submission

This is an update to `sfa`, from the current CRAN version 1.1.5 to 1.2.0.

There are **no reverse dependencies on CRAN**, so no other package is affected
by anything below.

### On the short interval since 1.1.5

1.1.5 was published on 2026-08-23, so `R CMD check --as-cran` raises the
"Days since last update" note, and I am aware that the policy asks for updates
to be spaced out. I am submitting this soon anyway because **1.1.5 silently
returns wrong results for one of its models**, and I would rather that version
were on CRAN for days than for months.

`sfm(model_name = "NE")` and `"NGE"` in 1.1.5 can return a **positive**
log-likelihood with `sigma_u` driven to zero, as an ordinary fitted object with
no error and no warning. `log Phi(z)` and an exponential tilt both diverge like
`z^2/2` with opposite signs as `sigma_u -> 0`, and at the scales the optimiser
visits their sum is catastrophic cancellation returning rounding noise — which
the optimiser then maximises by running `sigma_u` to its bound. Two examples
from a scan of 150 samples: `sigma_v = 587.4`, `sigma_u = 0`,
`logLik = +468992`; and `sigma_v = 3.4e15`, `logLik = +7.9e30`. Across a
12-cell design at 1,500 replications the rate was 0.74%, reaching 4.4% at
`lambda = 0.5`, `n = 100`.

1.2.0 does that cancellation analytically instead. The new expression agrees
with the old one to 8.5e-13 over 75,000 evaluations wherever the old one was
trustworthy, integrates to 1, and gave no failures in 2,200 fits.

The same release also fixes `"NGE"` aborting outright on about 7% of small
samples, and `nobs()` returning the number of rows supplied rather than the
number used, which made `BIC()` wrong for any fit on data with a missing value.

If you would prefer that I hold the new features and submit only the
corrections, or that I wait, I am happy to do either — please just say which.

### No breaking changes in this version

Unlike 1.1.5, nothing in this release changes the meaning of an existing
argument or `model_name`. Every addition below is a new argument defaulting to
the previous behaviour, or a new `model_name` value.

Two run-time warnings introduced in 1.1.5 -- for `psfm(model_name = "TFE")` and
`psfm(model_name = "GTRE")`, whose estimators changed in that version -- are
**retained** rather than removed, even though the one release cycle they were
promised for has now elapsed. 1.1.5 was the first release since 1.0.4 and is
only a few days old, so most users upgrading to 1.2.0 will be coming from
1.0.4 and meeting those changes for the first time.

### What is new

Four new model-fitting entry points, each for a model the package could not
previously express:

* `lcsfm()` -- the latent class stochastic frontier (Greene 2005; Orea and
  Kumbhakar 2004), with `n_class` coexisting technologies and posterior class
  probabilities.
* `selsfm()` -- Greene's (2010) frontier with a correction for sample
  selection, where the units in the sample are there for reasons correlated
  with the frontier's own noise. Two-step: probit, then maximum simulated
  likelihood. Its standard errors are conditional on the first stage and do not
  carry the Murphy-Topel correction, which `?selsfm` states plainly.
* `ivsfm()` -- endogenous regressors (Amsler, Prokhorov and Schmidt 2016,
  2017), with three estimators: full maximum likelihood, a two-step control
  function, and corrected 2SLS. To our knowledge no other CRAN package corrects
  for endogeneity in a stochastic frontier.
* `copsfm()` -- dependence between the noise and the inefficiency through a
  copula (Smith 2008), relaxing the independence assumption every other model
  in the package makes. `?copsfm` documents, with the numbers, that the
  dependence parameter is consistent but needs a large sample; the frontier
  slopes are unaffected.

A fifth specification, `lcsfm(model_name = "LCM_CN")`, the contaminated normal
frontier: every parameter is common across components except the noise scale,
so the noise has heavier tails without the frontier or the inefficiency varying.
Its composed density is a closed form -- a mixture of ordinary normal/half
normal densities sharing one `sigma_u` -- so it costs no extra integration.

**Model selection and specification testing**, which is the largest addition by
volume. The package offers fifteen cross-sectional inefficiency distributions
and previously gave the user nothing but AIC/BIC to choose among them:

* `TIC()` and `vuong()` -- Takeuchi's criterion and Vuong's non-nested test
  (Lai and Huang 2010). Neither assumes any candidate is correctly specified,
  which is the assumption in doubt when the choice is being made, and most of
  these distributions are not nested so the ordinary likelihood ratio test has
  no chi-square limit.
* `spec_test()` and `spec_test_all()` -- Papadopoulos and Parmeter's (2023)
  test of a distributional PAIR, computed from OLS residuals before any
  frontier is fitted. Twelve noise/inefficiency combinations.
* `lcsfm_homogeneity()` -- a test of a latent class fit against homogeneity.
* `sfma()` -- model averaging over inefficiency distributions (Parmeter, Wan
  and Zhang 2019), for when the honest answer is that the data does not
  identify one.

Each of these last three defaults to a bootstrap null rather than the published
asymptotic one, and `?spec_test`, `?lcsfm_homogeneity` and `?sfma` give the
measured size distortions that led to that choice. In two cases the asymptotic
null was badly mis-sized for the models this package fits -- 63.5 and 18.7
percent rejection at a nominal 5 percent -- because the published limits are
stated for restricted specifications the package does not impose.

New extractor and diagnostic functions: `efficiency()` (Battese-Coelli, JLMS
and modal predictors, on either scale of the dependent variable),
`meanefficiency()` (the model-implied `E[exp(-U)]` and supra-percentile means,
in closed form for seven distributions), `simulation_se()`, which reports how
much of a simulated-ML standard error is simulation noise rather than sampling
noise, and `pcomposed()`/`dcomposed()`, the distribution and density of the
composed error (Amsler, Schmidt and Tsay 2019), accurate in the far tail where
a copula argument needs it.

Further arguments: `psfm(mundlak = )` adds Mundlak adjustment terms for firm
effects correlated with the regressors (Karagiannis and Kellermann 2019), and
`lcsfm(penalty_c = )` maximises the modified likelihood the homogeneity test
requires.

Further arguments to `sfm()`: `weights`/`wscale`, `start_from` (seed a hard
model from a simpler fitted one, matched by parameter name), `scaling` for the
Wang and Schmidt (2002) scaling-property model, and `shapehet`, which is
labelled experimental in `?sfm` because its coefficients do not recover well in
testing. `vcov(type = "bhhh")` gives an outer-product-of-gradients covariance
that is defined where the Hessian is not, and `extract()` methods let
\pkg{texreg} render these fits into tables.

Six new panel estimators, all classical (non-maximum-likelihood) and none
assuming a distribution for the inefficiency term:

* `model_name = "CSS"` -- Cornwell, Schmidt and Sickles (1990)
* `model_name = "LS"` -- Lee and Schmidt (1993)
* `model_name = "KSS"` -- Kneip, Sickles and Song (2012)
* `model_name = "SSRE"` / `"SSCRE"` -- the random-effects and
  correlated-random-effects members of the Schmidt and Sickles (1984) family,
  whose within estimator the package already had as `"SSFE"`

Like the existing `"SSFE"`, none of these is maximum likelihood, so they carry
no optimisation object and `logLik()`/`AIC()`/`BIC()` return `NA` with a
warning rather than a number that would not mean what it appears to.

Heteroskedasticity in more than one error component, as named formula
arguments to `sfm()`:

* `vhet` for the noise scale, `uhet` for the inefficiency scale, `muhet` for
  the pre-truncation mean -- the last being Battese and Coelli (1995).

Also `z_link` for `sfm()` and `ttsfm()`, `marginal_effects()` for panel fits,
and `efficiency_ci()` for Horrace and Schmidt (1996) intervals.

### Changes that alter numeric output

Five. Three are corrections to defects; two change the simulation draws used
by the panel estimators, and are listed here because they move results for
existing users even though no argument changed meaning.

* `sfm(model_name = "NE")` and `"NGE"` could return a **positive**
  log-likelihood with `sigma_u` driven to zero. `log Phi(z)` and an
  exponential tilt both diverge like `z^2/2` with opposite signs as
  `sigma_u -> 0`, and their sum was catastrophic cancellation returning
  rounding noise, which the optimiser then maximised. Measured across a
  12-cell design at 1,500 replications the rate was 0.74%, reaching 4.4% in
  the worst cell. The cancellation is now done analytically.

* `nobs()` counted rows **supplied** rather than rows **used**, because it
  re-evaluated the `data` argument of the recorded call while the fitting code
  drops incomplete cases first. `BIC()` was therefore computed against the
  wrong `n` for any fit on data containing a missing value.

* `sfm(model_name = "NGE")` aborted outright on about 7% of small samples with
  "non-finite finite-difference value", because an out-of-domain guard
  returned `.Machine$double.xmax` and `optim()` differences the objective for
  its gradient.

* `psfm()` gave **every firm the same simulation draws**. One `R x 2` Halton
  block was built once and recycled across all `N` firms, in three separate
  places, so the negative correlation across observations that is half of
  Halton's advantage (Train 2002, p. 228) was absent, and every firm's
  simulation error was the same realization and could not average out as `N`
  grew. Each firm now gets its own contiguous block. Measured on 87 paired
  replications at N = 50/100/200 with identical seeds, `sigh` improves on 62
  of 87 (Wilcoxon p = 0.015) and `sigr` on 48 of 87 (p = 0.048), while the
  frontier slope -- the control, since it is not what the draws integrate
  over -- is unaffected at 46 of 87 (p = 0.41). Any simulated-ML panel fit
  changes slightly as a result.

* `rand.gtre` randomized the draws by a method that removed their point. It
  drew 9999 random permutations of Halton dimension 1 and kept whichever
  correlated least with dimension 2; permuting one column preserves each
  margin but randomizes the pairing, leaving joint coverage no better than
  random. Measured joint discrepancy at `R = 150`: 0.0217 before the shuffle,
  0.0550 after, against 0.1400 for purely random pairing. It now applies a
  uniform shift modulo 1 (Tuffin 1996; Train 9.3.4), which moves the lattice
  without disturbing its structure -- 0.0250 on the same measurement. Results
  change for anyone who passed `rand.gtre`.

### Dependencies

Unchanged. No package was added to `Imports` or `Suggests` in this version,
and none was removed.

## Vignette build time

Recorded here because it was the blocker on the 1.1.5 submission. The single
vignette is unchanged in this release and still renders in **13.9 seconds**
locally. The six new estimators are documented in `?psfm` with runnable
examples rather than in the vignette, deliberately, to keep the build time
where Uwe Ligges asked for it.

## Test environments

* local macOS 26.5 (aarch64-apple-darwin20), R 4.5.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (R-devel and R-release)

## R CMD check results

`R CMD check --as-cran --run-donttest`, R 4.5.2 on macOS 26.5:
**0 errors | 0 warnings | 3 notes.**

1. `checking CRAN incoming feasibility ... NOTE`, reporting the days since
   1.1.5. This is the point addressed at the top of this file.

2. `checking examples ... NOTE`, listing five examples over the 5-second
   guideline: `simulation_se` (~16 s), `zsfm` (9.9 s), `PL80_MVTN` (7.8 s),
   `sfa_diagnostics` (5.5 s) and `psfm` (5.0 s). `simulation_se()` estimates
   how much of a standard error is simulation noise by REFITTING the model K
   times with independent randomizations, so its example necessarily costs
   K + 1 simulated-ML fits; it was cut from 46 s by halving the sample and
   using the smallest K that still shows a spread. All four are inside `\donttest{}`; they fit models by
   simulated maximum likelihood over Halton draws, which is what costs the
   time. The `psfm` example was reduced from 20.4 s for this release, by
   shrinking its panel to 60 firms over 6 periods and passing
   `halton_num = 50` -- the same treatment the vignette received for 1.1.5.
   It is close to the floor: the `TRE_Z` fit it demonstrates becomes unstable
   below about 60 firms.

3. `checking HTML version of manual ... NOTE`, reporting that HTML Tidy is not
   recent enough and that package `V8` is unavailable, so those two sub-checks
   are skipped rather than failed. This is a property of the check machine
   rather than the package.

A fourth note, `checking for future file timestamps ... NOTE` / `unable to
verify current time`, appears intermittently here when the clock-check web
service is unreachable. It is also a property of the machine.

## Notes for the reviewer

* `NEWS.md` is long for one version. This release adds six estimators and
  several arguments, and the entries record the measurements behind the
  numerical fixes rather than only naming them.

* Tests needing a statistically meaningful sample size are behind
  `skip_on_cran()`, so the suite run during CRAN's checks is limited to fast
  structural tests. One further model (`ttsfm(model_name = "TTHN")`) is
  skipped unless `SFA_TEST_SLOW` is set, because it is much slower than the
  other estimators.

* Several examples are wrapped in `\donttest{}` because they fit models by
  simulated maximum likelihood over Halton draws, or by kernel regression with
  bandwidth cross-validation, and take longer than the 5-second guideline.
  They are checked with `--run-donttest` before submission, and the whole
  examples stage runs in 46 s.

* Every example using `PSopt = TRUE` now passes `rand.psoptim`. The
  particle-swarm stage draws from the session RNG, so without a seed the
  printed results change from build to build.

* `model_name = "KSS"` requires a balanced panel, which is what the estimator
  is defined on; it stops with a message naming `"CSS"` and `"LS"` as the
  unbalanced-panel alternatives rather than silently fitting something else.
