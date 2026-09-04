## Submission

This is an update to `sfa`, from the current CRAN version 1.1.5 to 1.2.0.

There are **no reverse dependencies on CRAN**, so no other package is affected
by anything below.

### Why now, so soon after 1.1.5

1.1.5 was published on 2026-08-23. `R CMD check --as-cran` no longer raises the
"Days since last update" note, but the interval is still short by the spirit of
the policy and I want to be explicit about the reason.

**1.1.5 can silently return wrong results for two of its models.**
`sfm(model_name = "NE")` and `"NGE"` can return a **positive** log-likelihood
with `sigma_u` driven to zero, as an ordinary fitted object with no error and
no warning. `log Phi(z)` and an exponential tilt both diverge like `z^2/2` with
opposite signs as `sigma_u -> 0`, and at the scales the optimiser visits their
sum is a catastrophic cancellation returning rounding noise -- which the
optimiser then maximises by running `sigma_u` to its bound. Across a 12-cell
design at 1,500 replications the rate was 0.74%, reaching 4.4% in the worst
cell. 1.2.0 performs that cancellation analytically instead.

Repairing it at one boundary opened the mirror image at the other, in this
release's own development rather than in anything published: the closed form
returns `log Phi(z) + z^2/2`, and for `eps < 0` with `sigma_v -> 0` the caller
then subtracted a term of the same magnitude, ~5e15, where consecutive doubles
are about 1 apart. That is fixed here too, by doing the second subtraction
analytically as well. Both were found by the same root-n convergence study, and
`"NE"` now passes it on the mean basis at 1000 replications per sample size
(slopes -1.02 to -1.04 against a target of -1, every R-squared at least 0.989,
no failures in 5000 fits) where before it did not. `NEWS.md` has the
derivations.

Two further corrections in the same release: `"NGE"` aborted outright on about
7% of small samples, and `nobs()` returned rows **supplied** rather than rows
**used**, which made `BIC()` wrong for any fit on data containing a missing
value.

I would rather 1.1.5 were the current version for weeks than for months, which
is why the corrections and the new features are submitted together rather than
the corrections alone.

### No intentional API-breaking changes in this version

Nothing in this release changes the meaning of an existing argument or
`model_name`. Every addition is a new argument defaulting to the previous
behaviour, or a new `model_name` value.

Existing calls can nonetheless return **different numbers**, and I would rather
flag that than have it noticed. Two changes to the simulation draws used by the
panel estimators move results for existing users:

* `psfm()` gave every firm the **same** simulation draws -- one Halton block
  was built once and recycled across all firms, so the negative correlation
  across observations that is half of Halton's advantage was absent and each
  firm's simulation error was the same realisation. Each firm now gets its own
  block. Any simulated-ML panel fit changes slightly.
* `rand.gtre` randomised the draws by a method that removed their point,
  permuting one Halton column and so randomising the pairing. It now applies a
  uniform shift modulo 1 (Tuffin 1996), which moves the lattice without
  disturbing its structure. Results change for anyone who passed `rand.gtre`.

`NEWS.md` gives the measurements behind both.

Two run-time warnings introduced in 1.1.5 -- for `psfm(model_name = "TFE")` and
`psfm(model_name = "GTRE")`, whose estimators changed in that version -- are
**retained** rather than removed, even though the one release cycle they were
promised for has elapsed. 1.1.5 was the first release since 1.0.4 and is only
two weeks old, so most users upgrading to 1.2.0 will be coming from 1.0.4 and
meeting those changes for the first time.

### What is new

Summarised here; `NEWS.md` has the detail.

**Four new model-fitting entry points**, each for a model the package could not
previously express: `lcsfm()` (latent class, Greene 2005; Orea and Kumbhakar
2004), `selsfm()` (sample selection, Greene 2010), `ivsfm()` (endogenous
regressors, Amsler, Prokhorov and Schmidt 2016) and `copsfm()` (copula
dependence between the two error components).

**Five new panel estimators** in `psfm()`, all classical rather than maximum
likelihood, and none assuming a distribution for the inefficiency term:
`"CSS"` (Cornwell, Schmidt and Sickles 1990), `"LS"` (Lee and Schmidt 1993),
`"KSS"` (Kneip, Sickles and Song 2012), and `"SSRE"` / `"SSCRE"`, the
random-effects and correlated-random-effects members of the Schmidt and Sickles
(1984) family whose within estimator the package already had as `"SSFE"`. Like
`"SSFE"`, none is maximum likelihood, so they carry no optimisation object and
`logLik()`/`AIC()`/`BIC()` return `NA` with a warning rather than a number that
would not mean what it appears to.

**Model selection and specification testing**, the largest addition by volume.
The package offers fifteen cross-sectional inefficiency distributions and
previously gave the user nothing but AIC/BIC to choose among them: `TIC()` and
`vuong()`, `spec_test()`/`spec_test_all()`, `lcsfm_homogeneity()`, and `sfma()`
for averaging over distributions rather than selecting one. The last three
default to a **bootstrap** null rather than the published asymptotic one, and
`?spec_test`, `?lcsfm_homogeneity` and `?sfma` give the measured size
distortions that led to that choice -- in two cases the asymptotic null was
badly mis-sized for the models this package fits, because the published limits
are stated for restricted specifications the package does not impose.

**A second sample-selection model.** `selsfm(model_name = "kts")` fits Kumbhakar,
Tsionas and Sipilainen (2009): two technologies, each with its own frontier and
scales, where the technology choice depends on inefficiency itself rather than
on the noise. Estimated by single-step maximum likelihood, because neither
two-step order is available -- the choice equation cannot be a probit when the
inefficiency entering it is unobserved.

**More copula families.** `copsfm()` offers fifteen where it offered two: Frank,
Clayton, Gumbel and Joe with rotations, alongside Gaussian and FGM. Each density
is checked against the second mixed partial of its own CDF. Most of them do not
recover their dependence parameter -- 36 to 60 percent of fits return the
independence boundary on data generated from that same family -- and `copsfm()`
warns, quoting the measured rate, rather than leaving a user to infer it from an
implausible estimate. `?copsfm` gives the table.

**Robustness diagnostics** for the divergence estimators `sfm()` already
offered: `hscore()`/`hscore_select()` and `calibrate_c()` to choose the tuning
parameter rather than supply it by hand, `density_weights()` for the weight each
observation receives, and `influence_sfa()` for the influence function of a fit.

**Other additions**: heteroskedasticity in more than one error component via
`vhet`/`uhet`/`muhet`; `z_link` for `sfm()` and `ttsfm()`; `efficiency()`,
`meanefficiency()`, `efficiency_ci()`, `marginal_effects()`, `simulation_se()`,
`pcomposed()`/`dcomposed()`; `weights`, `start_from`, `scaling` and an
experimental `shapehet` on `sfm()`; `vcov(type = "bhhh")`; and `extract()`
methods for \pkg{texreg}.

### Dependencies

Unchanged. No package was added to `Imports` or `Suggests` in this version, and
none was removed.

## Vignette build time

Recorded here because it was the blocker on the 1.1.5 submission. The single
vignette is unchanged in this release and still renders in **13.9 seconds**
locally. The five new panel estimators are documented in `?psfm` with runnable
examples rather than in the vignette, deliberately, to keep the build time where
Uwe Ligges asked for it.

## Test environments

* local macOS 26.5 (aarch64-apple-darwin20), R 4.5.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (R-devel and R-release)

## R CMD check results

`R CMD check --as-cran --run-donttest`, R 4.5.2 on macOS 26.5:
**0 errors | 0 warnings | 2 notes**, both of which are properties of the check
machine rather than of the package.

1. `checking CRAN incoming feasibility ... NOTE` --- reports the maintainer
   address, and that the package has a `VignetteBuilder` field but no prebuilt
   vignette index. The second half is an artefact of checking with
   `--no-build-vignettes`: pandoc is not installed on this machine, so the
   vignette cannot be built here. It builds on all five GitHub Actions
   platforms and on win-builder.

   The "Days since last update" note is **not** raised. It was when this file
   was first written; 1.1.5 is now twelve days old.

2. `checking top-level files ... NOTE` --- `README.md` and `NEWS.md` cannot be
   checked without pandoc. Same cause as above.

Two `WARNING`s also appear locally, both saying that `inst/doc` does not exist
and that `intro_to_psfm.Rmd` has no rendered output. Both are the same missing
pandoc. They do not appear where the vignette can be built.

`checking examples ... OK` and `checking examples with --run-donttest ... OK`.
The examples that dominate that stage, measured from `sfa-Ex.timings`:

\tabular{ll}{
 `influence_sfa` \tab 76 s \cr
 `simulation_se` \tab 38 s \cr
 `zsfm` \tab 23 s \cr
 `PL80_MVTN` \tab 15 s \cr
 `lcsfm_homogeneity` \tab 14 s
}

All are inside `\donttest{}`. They fit models by simulated maximum likelihood
over Halton draws, by quadrature per observation, or by bootstrap, which is
what costs the time. `simulation_se()` estimates how much of a standard error
is simulation noise by REFITTING the model K times with independent
randomizations, so its example necessarily costs K + 1 simulated-ML fits; it
was cut from 46 s by halving the sample and using the smallest K that still
shows a spread. `influence_sfa`'s example fits a Student-t noise model, which
is multi-start because that likelihood is bimodal in its degrees of freedom;
below n = 120 its Hessian stops inverting, so that is the floor rather than a
choice.

Two examples were cut for this submission after measuring rather than
estimating: `copsfm` from **239 s to 7 s**, by moving its example from n = 4000
at the default 128 quadrature nodes to n = 600 at 64 -- the release notes
already record that the quadrature is converged by 64, so nothing is lost --
and `influence_sfa` from 116 s to 76 s. The whole `--run-donttest` stage falls
from 511 s to about 240 s as a result.

Not reproduced locally, and expected on the submission machine:
`checking HTML version of manual ... NOTE`, reporting that HTML Tidy is not
recent enough and that package `V8` is unavailable, so those two sub-checks are
skipped rather than failed; and `checking for future file timestamps ... NOTE`
when the clock-check web service is unreachable. Both are properties of the
machine. The check above was run with `--no-manual`, so the first did not
arise.

## Notes for the reviewer

* `NEWS.md` is long for one version, and deliberately so: this release adds
  four entry points, five panel estimators and a number of arguments, and the
  entries record the measurements behind the numerical fixes rather than only
  naming them. Anything cut from this letter for length is there.

* Tests needing a statistically meaningful sample size are behind
  `skip_on_cran()`, so the suite run during CRAN's checks is limited to fast
  structural tests. One further model (`ttsfm(model_name = "TTHN")`) is
  skipped unless `SFA_TEST_SLOW` is set, because it is much slower than the
  other estimators.

* Several examples are wrapped in `\donttest{}` because they fit models by
  simulated maximum likelihood over Halton draws, or by kernel regression with
  bandwidth cross-validation, and take longer than the 5-second guideline.
  They are checked with `--run-donttest` before submission, which takes 292 s.

* Every example using `PSopt = TRUE` now passes `rand.psoptim`. The
  particle-swarm stage draws from the session RNG, so without a seed the
  printed results change from build to build.

* `model_name = "KSS"` requires a balanced panel, which is what the estimator
  is defined on; it stops with a message naming `"CSS"` and `"LS"` as the
  unbalanced-panel alternatives rather than silently fitting something else.
