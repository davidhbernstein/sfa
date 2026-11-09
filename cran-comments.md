## Submission

An update to `sfa`, from CRAN's current 1.2.0 to 1.3.0. There are **no reverse
dependencies on CRAN**, so no other package is affected.

### What the release is

1.3.0 adds estimators and fixes defects in released code. Most of those defects
returned a wrong or irreproducible result without an error or a warning, and
several fixes change estimates. `NEWS.md` groups them so that a user can see
what to re-run:

* **Breaking changes** (six items). For example, `ttsfm()` now reports its
  scales on their natural scale under their own names, and `Method` accepts
  only `"L-BFGS-B"`, because other methods dropped the bounds that keep the
  scales positive.
* **Changes that alter results.** Among them: `copsfm(inefdec = FALSE)` maximised
  the density of `v - u` instead of `v + u`; `sfm(estimator = "cols")` returned
  cost-frontier slopes with the wrong sign; a transformed response such as
  `log(y) ~ ...` was fitted untransformed; `"NW"` and `"NLN"` are now fitted by
  quadrature instead of simulation; and several fits that stopped short of a
  maximum, or whose result depended on the platform, now reach it.
* **New features**, including `sfm(estimator = "acols")` and `"cmle"`,
  `skewness_decomp()`, the `"NGB2"` and `"NB"` specifications, noise
  distributions for `copsfm()`, and sandwich and clustered covariances from
  `vcov()`.
* **Bug fixes** to standard errors, efficiency prediction, failures on valid
  input, and messages.

### Dependency changes since 1.2.0

`pracma`, `MASS` and `tmvtnorm` are no longer used. `tmvtnorm` is replaced by
`mvtnorm`, on which it was built; `tmvtnorm::ptmvnorm()` over an orthant is
exactly `mvtnorm::pmvnorm()`, and the efficiency scores are bitwise identical.
`erfinv()` and `ginv()` are internal copies of the same computation, verified
bitwise identical.

`gsl` and `pbapply` move from `Imports` to `Suggests`, so installing from
source no longer needs the GNU Scientific Library. `gsl` supplied only the
parabolic cylinder function for `"NG"`/`"NNAK"`, which is now computed in R
and agrees with it to about 1e-10; the tests keep `gsl` as an independent
reference. `pbapply` provided only an optional progress bar, already guarded by
`requireNamespace()`.

`lpSolve` leaves `Suggests` and `DEA` enters it. `npsfm(method = "SZ")` now
runs its DEA step through `DEA::dea()`, which the maintainer of this package
also maintains. The code path is guarded by `requireNamespace()`, and the two
implementations agreed to 6e-12 across every returns-to-scale setting before
the swap.

## R CMD check results

<!-- FILL ON THE RELEASE COMMIT (checklist item 10). Nothing below this line
     until "Test environments" has been measured on 1.3.0 yet. -->

One command, run once, on the tarball being submitted:

```
R CMD check --as-cran sfa_1.3.0.tar.gz
```

with `NOT_CRAN=true` set, R 4.5.2 on macOS (aarch64), on a tarball built
**with** the vignette. Result: **TO FILL**. The expected result is
0 errors | 0 warnings | 1 note, the note being a property of the check machine
rather than of the package:

1. `checking HTML version of manual ... NOTE`: HTML Tidy on this machine is
   not recent enough and package `V8` is unavailable, so the HTML-validation
   and math-rendering sub-checks are skipped rather than failed.

**A note on the command, because it changes the result.** `--as-cran` already
runs the `\donttest{}` examples, as a separate `checking examples with
--run-donttest` stage. Passing `--run-donttest` explicitly *in addition* folds
them into `checking examples` instead, which then exceeds the 5-second
guideline and raises a second note. That note is an artefact of the redundant
flag, measured on 1.2.1; the command above is the one reported.

## Check time

`checking tests ... TO FILL` under `NOT_CRAN=true`, which runs the Monte Carlo
and bootstrap validations that CRAN skips. Under CRAN's own conditions the
stage is much shorter, because the tests that need a statistically meaningful
sample size are behind `skip_on_cran()`. The warnings in that run are
deliberate diagnostics exercised by the tests that exist to fire them
(boundary reports, wrong-skew reports, rename notices), together with warnings
raised by `plm` and by `optim()`'s numerical Hessian. None accompanies a failed
expectation.

`checking examples` and `checking examples with --run-donttest`: TO FILL.

## Test environments

* local: R 4.5.2, macOS (aarch64)
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macos-latest,
  windows-latest
* win-builder (R-devel and R-release)

## Notes for the reviewer

* Examples are wrapped in `\donttest{}` where they fit models by simulated
  maximum likelihood over Halton draws, by quadrature, or by kernel regression
  with bandwidth cross-validation, and so exceed five seconds. They are checked
  with `--run-donttest` before submission, as above.

* Tests needing a statistically meaningful sample size are behind
  `skip_on_cran()`. One model, `ttsfm(model_name = "TTHN")`, is skipped unless
  `SFA_TEST_SLOW` is set, because a single fit can take an hour.

* Every example using `PSopt = TRUE` passes `rand.psoptim`. The particle-swarm
  stage draws from the session RNG, so without a seed the printed results
  change between builds.

* `model_name = "KSS"` requires a balanced panel, which is what the estimator
  is defined on. It stops with a message naming `"CSS"` and `"LS"` as the
  unbalanced-panel alternatives rather than silently fitting something else.
