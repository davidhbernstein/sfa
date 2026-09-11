## Submission

An update to `sfa`, from CRAN's current 1.2.0 to 1.2.1. There are **no reverse
dependencies on CRAN**, so no other package is affected.

### Why so soon after 1.2.0

1.2.0 was published three days ago, and `R CMD check --as-cran` raises the
"Days since last update" note accordingly. The reason is a defect in 1.2.0
that returns wrong estimates without any error or warning.

`copsfm(inefdec = FALSE)` -- the cost-frontier orientation of the new copula
entry point -- composed the density with the sign applied twice, so it
maximised the density of `v - u` where it needed the density of `v + u`. The
fit does not fail; it returns an ordinary-looking result with the frontier
slopes intact and the scales wrong. On clean cost data with `sigma_u = 1` and
`sigma_v = 0.4` it reported `sigma_u = 0.031` and `sigma_v = 0.69`, with the
intercept 0.8 too high -- which reads as "no inefficiency found".

Production fits, the default, were never affected, and are what the existing
tests covered. The repair is one line; the regression test that pins it checks
the composed density against its closed form in **both** orientations.

Two further defects in already-released code are fixed in the same submission,
both silent:

* A parameter converging **onto a bound** cost `copsfm()` every standard error
  in the fit. The likelihood refused out-of-range draws with
  `.Machine$double.xmax`; `optim()` differences the objective to form its
  gradient, and differencing 1.8e308 overflows, so the final stage aborted and
  all standard errors came back `NA`. Now a large finite penalty.
* `sfm(estimator = "cols")` used the normal/half-normal efficiency posterior
  regardless of `model_name`, so `"NE"` and `"NG"` fits received the wrong
  `exp_u_hat`. Parameter estimates were unaffected.

The release also adds three inefficiency/noise specifications and rewrites the
package Description. Those are in `NEWS.md`; they are not the reason for the
timing.

## R CMD check results

One command, run once, on the tarball being submitted:

```
R CMD check --as-cran sfa_1.2.1.tar.gz
```

with `NOT_CRAN=true` set, R 4.5.2 on macOS 26.5, on a tarball built **with**
the vignette. Result: **0 errors | 0 warnings | 2 notes**.

1. `checking CRAN incoming feasibility ... NOTE` -- "Days since last update",
   the interval explained above.

2. `checking HTML version of manual ... NOTE` -- HTML Tidy on this machine is
   not recent enough and package `V8` is unavailable, so the HTML-validation
   and math-rendering sub-checks are skipped rather than failed. A property of
   the check machine, not of the package.

Every other stage is `OK`, including `checking examples`, `checking examples
with --run-donttest`, `checking tests`, `checking top-level files`, `checking
files in 'vignettes'`, `checking package vignettes`, `checking re-building of
vignette outputs` and `checking PDF version of manual`.

**A note on the command, because it changes the result.** `--as-cran` already
runs the `\donttest{}` examples, as a separate `checking examples with
--run-donttest` stage. Passing `--run-donttest` explicitly *in addition* folds
them into `checking examples` instead, which then exceeds the 5-second
guideline and raises a third note listing thirteen examples. That note is an
artefact of the redundant flag: measured on this tarball, `--as-cran` alone
gives `checking examples ... [15s/16s] OK` followed by `checking examples with
--run-donttest ... [179s/176s] OK`, while `--as-cran --run-donttest` gives
`checking examples ... [200s/197s] NOTE`. Same examples, same machine, same
tarball. The command above is the one reported.

## Check time

`checking tests ... [33m/33m] OK` under `NOT_CRAN=true`, which runs the Monte
Carlo and bootstrap validations that CRAN skips: `FAIL 0 | WARN 16 | SKIP 8 |
PASS 3694`. Under CRAN's own conditions that stage is about a minute, because
the tests needing a statistically meaningful sample size are behind
`skip_on_cran()`. The 16 warnings are deliberate diagnostics being exercised by
the tests that exist to fire them -- boundary reports from Greene's true fixed
effects likelihood and from `GTRE`, the wrong-skew report from `npsfm("FLW")`,
and the `model_name = "TFE"` rename notice -- together with warnings raised by
`plm` and by `optim()`'s numerical Hessian stepping outside its own box. None
accompanies a failed expectation.

`checking examples` is 16 seconds and `checking examples with --run-donttest`
is 176 seconds. The examples that dominate the second are `influence_sfa`
(54.9s), `simulation_se` (28.5s) and `zsfm` (19.7s); all are inside
`\donttest{}` because they fit models by simulated maximum likelihood over
Halton draws, by quadrature, or by kernel regression with bandwidth
cross-validation.

## Test environments

* local: R 4.5.2, macOS 26.5 (aarch64)
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
