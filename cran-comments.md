## Submission

This is an update to `sfa`, from the current CRAN version 1.1.5 to 1.2.0.

There are **no reverse dependencies on CRAN**, so no other package is affected
by anything below.

### On the short interval since 1.1.5

`R CMD check --as-cran` notes "Days since last update: 5", and I am aware that
the policy asks for updates to be spaced out. I am submitting this soon anyway
because **1.1.5 silently returns wrong results for one of its models**, and I
would rather that version were on CRAN for days than for months.

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

Three, all of them corrections:

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

* local macOS 15 (aarch64-apple-darwin20), R 4.5.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (R-devel and R-release)

## R CMD check results

`R CMD check --as-cran --run-donttest`, R 4.5.2 on macOS 15:
**0 errors | 0 warnings | 3 notes.**

1. `checking CRAN incoming feasibility ... NOTE`, reporting "Days since last
   update: 5". This is the point addressed at the top of this file.

2. `checking examples ... NOTE`, listing four examples over the 5-second
   guideline: `zsfm` (10.4 s), `PL80_MVTN` (7.6 s), `sfa_diagnostics` (5.5 s)
   and `psfm` (5.2 s). All four are inside `\donttest{}`; they fit models by
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
