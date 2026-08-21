## Submission

This is a feature update to `sfa`, from the current CRAN version 1.0.4 to
1.1.3. It bundles the work of several development versions (1.1.0-1.1.3) that
were never submitted, so `NEWS.md` is correspondingly long.

Please note one deliberate **breaking change**, described first in `NEWS.md`:
`psfm(model_name = "TFE")` previously selected Chen, Schmidt and Wang's (2014)
within maximum-likelihood estimator and now selects Greene's (2005) true fixed
effects estimator, which is the standard meaning of that name in the
literature. The previous estimator is unchanged and available as
`model_name = "TFE_WMLE"`. Because existing scripts would otherwise silently
get a different estimator, `psfm()` emits a warning whenever `"TFE"` is used;
the warning is intended to remain for one release cycle. There are no reverse
dependencies on CRAN, so no other package is affected.

The package's dependency on `frontier` has been removed (along with eight
others), and the R requirement has been lowered from `R (>= 4.4.0)` to
`R (>= 4.0.0)`.

This version adds a fifth entry point, `npsfm()`, for nonparametric stochastic
frontier models. It needs kernel regression from `np`, and its `"SZ"` method
additionally needs `Benchmarking`. Both are declared in **Suggests** rather than
Imports, since nothing else in the package uses them: `npsfm()` tests for each
at run time and stops with an install instruction if it is absent, and its
examples and tests are guarded with `requireNamespace()` / `skip_if_not_installed()`
so they are skipped rather than failing where those packages are unavailable.
`npsfm()` returns an object of class `"npsfareg"`, not `"sfareg"` — a
kernel-estimated frontier has no parameter vector with standard errors, so the
`coef`/`vcov`/`logLik` methods do not apply to it.

## Test environments

* local macOS 15 (aarch64-apple-darwin), R 4.5.2
* [to be filled in before submission: win-builder release + devel, R-hub]

## R CMD check results

0 errors | 0 warnings | 1-2 notes

Both notes are properties of the local check machine rather than the package,
and are not expected on CRAN's systems:

* `checking HTML version of manual ... NOTE`
  `Skipping checking HTML validation: 'tidy' doesn't look like recent enough`
  `HTML Tidy` / `Skipping checking math rendering: package 'V8' unavailable`.
  The local HTML Tidy and V8 installations are missing or out of date, so
  those two sub-checks are skipped rather than failed.

* `checking for future file timestamps ... NOTE` / `unable to verify current
  time`, which appears intermittently when the clock-check web service is
  unreachable.

## Notes for the reviewer

* A `testthat` suite is new in this version. The tests that require a
  statistically meaningful sample size are behind `skip_on_cran()`, so the
  suite run during CRAN's checks is limited to fast structural tests. One
  further model (`ttsfm(model_name = "TTHN")`) is skipped unless the
  environment variable `SFA_TEST_SLOW` is set, because it is currently much
  slower than the other estimators.

* Several examples are wrapped in `\donttest{}` because they fit models by
  simulated maximum likelihood over Halton draws and take longer than the
  5-second guideline.
