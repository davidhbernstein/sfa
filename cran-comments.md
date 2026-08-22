## Resubmission

This is a resubmission. The previous submission of 1.1.4 checked with
`Status: 1 NOTE` on both pretest machines (Windows R-devel and Debian
R-devel), with no errors or warnings. The NOTE was `checking CRAN incoming
feasibility`, raising two items, both in `README.md`. Both are now fixed:

* **`URI: LICENSE.md`, "possibly invalid file URI".** This was a genuine
  broken link. `README.md` referred to `LICENSE.md` with a relative path,
  but `LICENSE.md` is listed in `.Rbuildignore` and so is not in the tarball,
  leaving the link dangling wherever the shipped README is rendered. It now
  points at the copy in the repository, which is where the file actually
  lives.

* **`https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html`, "Timeout was
  reached".** The URL is valid -- it returns 200 in about 0.15 s from here --
  and both machines reported a connection timeout rather than a bad response
  (21 s on Windows, 60 s on Debian). Rather than ask you to re-run it, the
  licence badge now links to CRAN's own copy of the GPL-2 text,
  <https://cran.r-project.org/web/licenses/GPL-2>, which cannot time out from
  a CRAN machine.

No other file changed, and there is no change to any R code, to `NAMESPACE`,
or to the documented interface.

## Submission

This is a feature update to `sfa`, from the current CRAN version 1.0.4 to
1.1.4. It bundles the work of five development versions (1.1.0-1.1.4) that
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
additionally needs `Benchmarking`. Both are declared in **Suggests** rather
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
