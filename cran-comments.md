## Resubmission

This is a resubmission, addressing Uwe Ligges's request to reduce the
vignette build time. Both pretests of 1.1.4 checked `Status: OK`, and
`checking re-building of vignette outputs` was reported at 10 minutes.

The package has one vignette, `intro_to_psfm`, which fits panel stochastic
frontier models by simulated maximum likelihood over Halton draws. It now
uses toy data and few iterations throughout, following your first two
suggestions:

* the simulated panel is 70 firms over 6 periods, down from 100 over 10;

* the simulated-ML fits pass `halton_num = 50`, where the default is
  `ceiling(sqrt(nrow(data))) + 100` -- about 130 draws at the old data size.
  The number of draws is the dominant cost in these likelihoods;

* the parametric-bootstrap example runs 5 replications on a 30-firm panel,
  down from 10 replications on 60 firms. Each replication refits the model,
  so this was the single most expensive chunk;

* the second simulated data set has been dropped. It was generated with
  arguments identical to the first, so the `GTRE_Z` section now reuses the
  panel already in hand.

Measured here, rendering the vignette went from 68.8 s to 13.2 s, a factor
of 5.2. I have not shrunk it further because the four-component GTRE model
becomes genuinely unstable on very short panels -- below about 70 firms its
efficiency-posterior step can fail to invert on some random seeds -- and I
would rather the vignette be fast than be fast and fragile. The chosen size
was checked against 12 seeds with no failures.

Two related points, in the interest of the vignette staying cheap and
predictable to check:

* every fit now sets `rand.gtre` and `rand.psoptim`. With `PSopt = TRUE` the
  particle-swarm stage draws from the session RNG, so the vignette's printed
  results previously changed from build to build. They are now reproducible.

* the vignette states explicitly that its settings are chosen to build
  quickly and that the estimates are not to be read as a serious fit.

No R code, `NAMESPACE`, or documented interface changed in this version.

## Submission

This is a feature update to `sfa`, from the current CRAN version 1.0.4 to
1.1.5. It bundles the work of six development versions (1.1.0-1.1.5) that
were never submitted, so `NEWS.md` is correspondingly long.

There are **no reverse dependencies on CRAN**, so no other package is affected
by anything below.

### Two deliberate breaking changes

Both are described at the top of their `NEWS.md` sections, and both warn at
run time rather than changing behaviour silently.

* `psfm(model_name = "TFE")` previously selected Chen, Schmidt and Wang's
  (2014) within maximum-likelihood estimator and now selects Greene's (2005)
  true fixed effects estimator, which is the standard meaning of that name in
  the literature. The previous estimator is unchanged and available as
  `model_name = "TFE_WMLE"`.

* `psfm(model_name = "GTRE")` now defaults to full information maximum
  likelihood rather than simulated maximum likelihood. The four ways of
  fitting the four-component GTRE model used to be four separate `model_name`
  values, which presented them as four models rather than four routes to one;
  they are now chosen with an `estimator` argument, in the same spirit as
  `sfm()`'s existing `estimator = c("mle", "cols")`. The old behaviour is
  `estimator = "sml"`.

In both cases existing scripts would otherwise silently get a different
estimator, so `psfm()` warns whenever the affected name is used without an
explicit choice. The warnings are intended to remain for one release cycle.

### New dependencies in Suggests

This version adds a fifth entry point, `npsfm()`, for nonparametric stochastic
frontier models. It needs kernel regression from `np`, and its `"SZ"` method
solves an output-oriented DEA as one linear program per unit using `lpSolve`.
Both are declared in **Suggests** rather
than **Imports**, since nothing else in the package uses them: `npsfm()` tests
for each at run time and stops with an install instruction if it is absent, and
its examples and tests are guarded with `requireNamespace()` /
`skip_if_not_installed()` so they are skipped rather than failing where those
packages are unavailable.

`npsfm()` returns an object of class `"npsfareg"`, not `"sfareg"` -- a
kernel-estimated frontier has no parameter vector with standard errors, so the
`coef`/`vcov`/`logLik` methods do not apply to it.

The dependency on `frontier` was removed in this cycle (along with eight
others), and the R requirement was lowered from `R (>= 4.4.0)` to
`R (>= 4.0.0)`.

## Test environments

* local macOS 15 (aarch64-apple-darwin20), R 4.5.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (R-devel and R-release)

## R CMD check results

0 errors | 0 warnings | 0-2 notes

Any notes seen locally are properties of the check machine rather than the
package, and are not expected on CRAN's systems:

* `checking HTML version of manual ... NOTE`, reporting that HTML Tidy is not
  recent enough and that package `V8` is unavailable, so those two sub-checks
  are skipped rather than failed.

* `checking for future file timestamps ... NOTE` / `unable to verify current
  time`, which appears intermittently when the clock-check web service is
  unreachable.

## Notes for the reviewer

* A `testthat` suite is new in this version. The tests that need a
  statistically meaningful sample size are behind `skip_on_cran()`, so the
  suite run during CRAN's checks is limited to fast structural tests. One
  further model (`ttsfm(model_name = "TTHN")`) is skipped unless the
  environment variable `SFA_TEST_SLOW` is set, because it is currently much
  slower than the other estimators.

* Several examples are wrapped in `\donttest{}` because they fit models by
  simulated maximum likelihood over Halton draws, or by kernel regression with
  bandwidth cross-validation, and take longer than the 5-second guideline.
  They are checked with `--run-donttest` before submission.
