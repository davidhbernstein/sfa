## Submission

An update to `sfa`, from CRAN's current 1.2.0 to 1.2.1. There are **no reverse
dependencies on CRAN**, so no other package is affected.

### Why so soon after 1.2.0

1.2.0 was published on 2026-09-07. The reason for an update within weeks is a
defect in 1.2.0 that returns wrong estimates without any error or warning.

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

Fifteen further defects in already-released code are fixed in the same submission,
all silent -- each returns a wrong or irreproducible result without an error:

* A parameter converging **onto a bound** cost `copsfm()` every standard error
  in the fit. The likelihood refused out-of-range draws with
  `.Machine$double.xmax`; `optim()` differences the objective to form its
  gradient, and differencing 1.8e308 overflows, so the final stage aborted and
  all standard errors came back `NA`. Now a large finite penalty.
* `sfm(estimator = "cols")` used the normal/half-normal efficiency posterior
  regardless of `model_name`, so `"NE"` and `"NG"` fits received the wrong
  `exp_u_hat`. Parameter estimates were unaffected.
* `sfm(estimator = "cols")` returned the frontier coefficients with the wrong
  **sign** on a cost frontier (`inefdec = FALSE`). The moment path fits on the
  production orientation, so a cost frontier is fitted as `-y = x'(-b) + ...`
  and `-b` was reported. Scale parameters were unaffected, which is why it was
  not caught. Cost-frontier `"cols"` fits made with 1.2.0 should be re-run.
* `psfm(model_name = "FD")`, `"TFE"` and `"TFE_WMLE"` were not reproducible,
  and each fit advanced the caller's random-number stream. An internal
  starting-value fit drew unseeded `runif()` starting values; two fits of
  identical data could differ by up to about 3e-4 in relative terms. The draws
  now use a fixed local seed and the caller's RNG state is restored.
* `psfm(model_name = "GTRE")` by simulated ML and `"GTRE_Z"` reported
  efficiency scores that differed between two fits of identical data, by up to
  0.4%, and advanced the caller's random-number stream. The scores come from
  `tmvtnorm::ptmvnorm()`, a randomized quasi-Monte Carlo integrator, which drew
  from the session RNG. It now runs under a fixed local seed and the caller's
  RNG state is restored. Parameter estimates were unaffected.
* `sfm(model_name = "NR")` computed the one-sided factor of its log-density as
  a difference of two nearly equal terms, which from z of about 7 overstated
  the log-density by up to about 7. It is now evaluated through the continued
  fraction for the Mills ratio, accurate to 1e-13; `exp_u_hat` likewise.
* Battese-Coelli efficiencies for `sfm()`'s `"NHN"`/`"NHN_Z"` and `psfm()`'s
  `"TRE"`/`"TRE_Z"` were a ratio of `1 - pnorm()` terms that underflowed for
  firms well above the frontier, reporting them at 0 instead of near 1. Now in
  logs; ordinary fits agree to 3e-15.
* `sfm(model_name = "NU")` floored a difference of normal CDFs that rounds to
  0 above the frontier, returning a log-density of -708 where the truth is
  -44; `"NHN"`/`"NHN_Z"` floored their density in levels. Both are now computed
  in logs, with estimates unchanged to within 8e-8 (NHN).
* `zsfm()`'s inefficient-regime log-density became `-Inf`, and the JLMS
  predictor in `zsfm()` and `lcsfm()` became `NaN`, for observations far from
  the frontier, through `log(dnorm())` and `dnorm()/pnorm()` in levels. Both now
  in logs; ordinary fits are unchanged.
* `zsfm(model_name = "ZISF_Z")` stopped in a spurious local optimum of its
  regime link in 7 of 50 simulated samples, 10 to 64 log-likelihood units below
  the best maximum and with the wrong regime slope. It now tries several link
  intercepts before optimizing, as `"ZISF"` already did; on the same samples
  no fit is more than 0.53 short.
* `ivsfm()` floored `pnorm()` at machine epsilon in its JLMS predictor, so a
  firm well above the frontier got a negative inefficiency and an efficiency
  above 1 (71 for a firm 20 above the frontier under `"C2SLS"`). Now in logs;
  estimates are unaffected.
* `ttsfm(model_name = "TTNE")` capped the exponents in its log-likelihood before
  forming it in levels. It is now a log-sum-exp, exact without a cap; fits are
  unchanged to about 1e-6 on simulated data.
* `ttsfm(model_name = "TTNE")`'s E[exp(-u) | eps], and so its M6 metric, used
  `0.5 * sig.v` where the closed form has `0.5 * sig.v^2`: 10.9% too high at
  `sigma_v = 0.3`. Corrected; checked against numerical integration.
* `sfm(model_name = "NR")`'s log-density formed a difference of two very large
  terms that lost its digits as `sigma_v` fell toward zero, reading too high,
  so fits ran to the `sigma_v` bound (reported -539.9 against an exact -585.1
  on one sample). Now formed exactly.
* `TIC()` returned -2.64e14 for a fit with a parameter on its bound, and
  `sfma(weights = "tic")` gave that model all the weight. It now refuses an
  indefinite Hessian or a non-positive penalty.

Eight further defects in released code made a call fail rather than return a
wrong answer:

* `psfm()` with `"PL80"`, `"BC92"`, `"K1990"`, `"K1990modified"` or `"SSFE"`
  stopped with "system is exactly singular" on a between-rank-deficient design
  such as `factor(year)` in an unbalanced panel. A random-effects starting
  regression none of these models uses was being run for them outside the
  collinearity guard; it is now run only for the models that read it.
* `psfm(keep_objective = TRUE)` stored a likelihood for `"GTRE"`, `"TRE"`,
  `"GTRE_Z"` and `"TRE_Z"`, but `influence_sfa()`, `vcov(type = "bhhh")`, the
  `sandwich` methods, `TIC()` and `vuong()` then failed on it. The retained
  likelihoods now return per-firm contributions.
* `psfm(model_name = "SSRE")` and `"SSCRE"` stopped with a bare LAPACK error
  on a singular random-effects fit. They now name the collinear columns, or,
  for `"SSCRE"` on an unbalanced panel, explain that this is a current
  limitation of that estimator and point to `"SSFE"`.
* `psfm(model_name = "GTRE_Z")` and `"TRE_Z"` stopped with "0 < ctrl$rhoend
  is not TRUE" when the random-effects regression seeding them estimated zero
  firm-effect variance: the starting step then began at zero scale, leaving
  `bobyqa()` no trust region. It now starts from a small positive scale.
* A starting value just outside its bound stopped a fit with "Starting values
  violate bounds" from `bobyqa()`; `psfm(model_name = "TRE")` on a boundary
  random-effects start fitted on macOS and failed on Linux and Windows. Only
  coordinates strictly outside their bounds are now moved inside; fits from
  valid starts are unchanged.
* `psfm(model_name = "GTRE_Z")` stopped when one firm's efficiency posterior
  could not be inverted, losing estimates that had converged. That firm now
  gets `NA` efficiencies with a warning.
* `psfm(model_name = "GTRE_Z")` stopped on Windows and Linux with "sigma must
  be positive definite" when a firm's posterior covariance lost positive
  definiteness to rounding. It is now repaired when the defect is
  rounding-sized; results where it was already valid are identical.
* `psfm(model_name = "SSCRE")` refused every unbalanced panel: the Mundlak
  means made `plm`'s Swamy-Arora step singular. It now uses the Swamy-Arora
  components of the model without the means, which on a balanced panel are
  exactly the ones it used before, so balanced results are unchanged.

The release also adds three inefficiency/noise specifications, a
wrong-skewness suite of two estimators and a diagnostic, and rewrites the
package Description. Those are in `NEWS.md`; they are not the reason for the
timing.

### Dependency changes

`lpSolve` leaves `Suggests` and `DEA` enters it. `npsfm(method = "SZ")` solved
its own output-oriented envelopment on `lpSolve`; it now calls `DEA::dea()`,
which the maintainer of this package also maintains. Both are `Suggests` and
the code path is guarded by `requireNamespace()`, so nothing in `sfa` requires
either package to be present. The two implementations agreed to 6e-12 across
every returns-to-scale setting before the swap.

`pracma` and `MASS` leave `Imports`, and `pbapply` moves from `Imports` to
`Suggests`. Each was used for one function: `erfinv()` and `ginv()` are now
internal copies of the same computation (results verified bitwise identical),
and `pbapply` only provided an optional progress bar that the code already
guarded with `requireNamespace()`.

## R CMD check results

One command, run once, on the tarball being submitted:

```
R CMD check --as-cran sfa_1.2.1.tar.gz
```

with `NOT_CRAN=true` set, R 4.5.2 on macOS 26.5.2 (aarch64), on a tarball built
**with** the vignette. Result: **0 errors | 0 warnings | 1 note**, a property
of the check machine rather than of the package:

1. `checking HTML version of manual ... NOTE` -- HTML Tidy on this machine is
   not recent enough and package `V8` is unavailable, so the HTML-validation
   and math-rendering sub-checks are skipped rather than failed.

An earlier run of this version also noted `checking for future file timestamps`
("unable to verify current time") when no time server could be reached; that
stage is `OK` in the run reported here.

Every other stage is `OK`, including `checking examples`, `checking examples
with --run-donttest`, `checking tests`, `checking top-level files`, `checking
files in 'vignettes'`, `checking package vignettes`, `checking re-building of
vignette outputs` and `checking PDF version of manual`.

**A note on the command, because it changes the result.** `--as-cran` already
runs the `\donttest{}` examples, as a separate `checking examples with
--run-donttest` stage. Passing `--run-donttest` explicitly *in addition* folds
them into `checking examples` instead, which then exceeds the 5-second
guideline and raises a second note listing seven examples. That note is an
artefact of the redundant flag. It was measured on an earlier build of this
version: there `--as-cran` alone gave `checking examples ... [11s/11s] OK`
followed by `checking examples with --run-donttest ... [127s/128s] OK`, while
`--as-cran --run-donttest` on the same tarball gave `checking examples ...
[126s/127s] NOTE`. On the tarball being submitted, `--as-cran` alone gives
`checking examples ... [11s/11s] OK` and `checking examples with --run-donttest
... [132s/133s] OK`. The command above is the one reported.

## Check time

`checking tests ... [25m/25m] OK` under `NOT_CRAN=true`, which runs the Monte
Carlo and bootstrap validations that CRAN skips: `FAIL 0 | WARN 16 | SKIP 9 |
PASS 3987`. Under CRAN's own conditions that stage is about a minute, because
the tests needing a statistically meaningful sample size are behind
`skip_on_cran()`. The 16 warnings are deliberate diagnostics being exercised by
the tests that exist to fire them -- boundary reports from Greene's true fixed
effects likelihood and from `GTRE`, the wrong-skew report from `npsfm("FLW")`,
and the `model_name = "TFE"` rename notice -- together with warnings raised by
`plm` and by `optim()`'s numerical Hessian stepping outside its own box. None
accompanies a failed expectation.

`checking examples` is 11 seconds and `checking examples with --run-donttest`
is 132 seconds. The examples that dominate the second are `influence_sfa`
(35.6s), `zsfm` (19.7s) and `simulation_se` (17.5s); all are inside
`\donttest{}` because they fit models by simulated maximum likelihood over
Halton draws, by quadrature, or by kernel regression with bandwidth
cross-validation.

## Test environments

* local: R 4.5.2, macOS 26.5.2 (aarch64)
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
