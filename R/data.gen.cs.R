## ------------------------------------------------------------------
## Cross-sectional data generator.
##
## See DATA_GENERATION_REFERENCE.md (project root) for the full mapping
## of every column here to the sfm()/zsfm()/ttsfm() model_name it's meant
## to test, the true parameter vector each model should recover, and an
## explicit list of which models are NOT yet covered by any column here.
## That file is the canonical reference; comments below are a shorter
## pointer to it, not a duplicate of it.
##
## `sig_w`: scale of the second one-sided component used only by the
## two-tier (TTNE/TTHN/TTNLS) columns below. Defaults to sig_u.
## ------------------------------------------------------------------
data_gen_cs <- function(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu, sig_w = sig_u){
if(!is.null(rand)){
set.seed(rand)}

n        <- N                                   ## Observations
x1       <- log(runif(n,2,10))                  ## x1 variable
x2       <- log(runif(n,3,50))                  ## x2 variable

## Normal half normal -- targets sfm(model_name="NHN") (and "NR", which is
## an alternative closed-form derivation of the same normal/half-normal
## composed error, so it reuses y_pcs too).
u        <- abs(rnorm(n,0,sig_u))               ## u_it half-normal error term
v        <- rnorm(    n,0,sig_v)

## Normal truncated normal -- targets sfm(model_name="NTN").
u_tn     <- rtruncnorm(n=n, mean = mu, sd=sig_u, a=0)

## NHN with Z -- targets sfm(model_name="NHN_Z", formula = y_pcs_z~x1+x2|z).
## sigma_u given z is exp(0.9 + 0.6*z) directly (NOT sqrt(exp(...)) -- see
## the reference doc for why this differs from psfm()'s GTRE_Z/TRE_Z
## convention). True z-coefficients to recover: intercept 0.9, slope 0.6.
z          <- runif(n,1,2)
uz         <- abs(rnorm(n, 0, exp(0.9 + 0.6*z)))       ## u_it half-normal error term with z var and cons

## Normal exponential -- targets sfm(model_name="NE").
phi      <- 1/sig_u
u_e      <- rexp(n,rate=phi)                    ## u_it exponential

## Normal exponential with Z -- targets sfm(model_name="NE_Z",
## formula = y_pcs_ez~x1+x2|z). Same exp(0.9+0.6*z) convention as uz above
## (NE_Z's sigma_u_fun in sfm.R uses the identical direct-exp link), just
## applied to an exponential rate instead of a half-normal sd.
uz_e     <- rexp(n, rate = 1/exp(0.9 + 0.6*z))

## T half T -- targets sfm(model_name="THT"). Shares the single "a"
## degrees-of-freedom parameter across the t-distributed noise and
## half-t inefficiency, matching THT's likelihood in sfm.R.
u_t      <- sig_u*abs(rt(n,a))                  ## u_it half-T error term
v_t      <- sig_v*    rt(n,a)

## Cauchy half Cauchy -- NOT currently matched to any sfm()/psfm() model
## (no "half-Cauchy" model_name exists). Kept for robustness/contaminated-
## noise sensitivity checks outside the main package's testable surface.
u_c      <- abs(rcauchy(n,0,sig_u))             ## u_it half-cauchy error term
v_c      <- rcauchy(    n,0,sig_v)

## Normal half uniform -- same status as the Cauchy block above: no
## matching model_name in the package, kept for exploratory/robustness use.
u_u      <- runif(n, min=0, max=sig_u)          ## u_it half-uniform error term
name     <- 1:N                                 ## name for individuals (cross-sectional: one obs each)
cons     <- cons                                ## constant

## ------------------------------------------------------------------
## Two-tier components (ttsfm()).
##
## Matches ttsfm.R's actual likelihoods: TTNE's fn() exponentiates
## alpha = e/sigu + sigv^2/(2*sigu^2) and
## a = sigv^2/(2*sigw^2) - e/sigw where e = Y - X*beta, which is exactly
## the closed-form composed-error density for
##   e = v + w - u,   v ~ N(0, sigv^2),  u ~ Exponential(mean = sigu),
##                                        w ~ Exponential(mean = sigw)
## (Polachek & Yoon 1987-style two-tier). TTHN's fn() (theta1/theta2/s/
## omega1/omega2 reparameterization) is the same e = v + w - u structure
## with u, w ~ Half-Normal(sigu), Half-Normal(sigw) instead. TTNLS
## (nonlinear least squares) uses E[u]=sigu, E[w]=sigw directly -- exactly
## the exponential case's first moments -- so it targets the same
## homoskedastic exponential columns as TTNE.
## ------------------------------------------------------------------

## Homoskedastic two-tier, exponential -- targets ttsfm(model_name="TTNE")
## and ttsfm(model_name="TTNLS"), formula = y_ttne~x1+x2 (no pipes).
w_tt     <- rexp(n, rate = 1/sig_w)             ## w_it exponential, second one-sided component
y_ttne   <-         cons + beta1*x1    + beta2*x2   + v + w_tt - u_e

## Homoskedastic two-tier, half-normal -- targets
## ttsfm(model_name="TTHN"), formula = y_tthn~x1+x2 (no pipes).
w_tt_hn  <- abs(rnorm(n, 0, sig_w))             ## w_it half-normal, second one-sided component
y_tthn   <-         cons + beta1*x1    + beta2*x2   + v + w_tt_hn - u

## Heteroskedastic (z/zp) two-tier, half-normal -- targets
## ttsfm(model_name="TTHN", formula = y_tthn_z~x1+x2|z|zp). sigma_u given z
## reuses the same exp(0.9+0.6*z) convention as uz/NHN_Z above; sigma_w
## given zp is exp(0.7+0.5*zp) (same direct-exp convention, both pipes are
## parameterized identically in ttsfm.R's TTHN fn). True coefficients to
## recover: sigma_u ~ z intercept 0.9 / slope 0.6, sigma_w ~ zp intercept
## 0.7 / slope 0.5.
zp        <- runif(n,1,2)
wz_hn     <- abs(rnorm(n, 0, exp(0.7 + 0.5*zp)))
y_tthn_z  <-         cons + beta1*x1    + beta2*x2   + v + wz_hn - uz

## ------------------------------------------------------------------
## Zero-inefficiency (zsfm()).
##
## ZISF's likelihood in zsfm.R is a two-component mixture,
## prob*N(0,sigv^2) + (1-prob)*[normal/half-normal composed-error density],
## i.e. each observation is drawn from an "efficient firm" regime (u=0,
## probability prob) or an "inefficient firm" regime (u ~ Half-Normal(sigu),
## probability 1-prob) -- the standard zero-inefficiency SFA formulation.
## ------------------------------------------------------------------

## ZISF (no z) -- targets zsfm(model_name="ZISF", formula=y_zisf~x1+x2).
## ZISF's own prob link is prob = exp(-abs(gamma)); using gamma0 = 0.3 here
## means the true recoverable gamma is 0.3 (up to the sign ambiguity
## |.| already bakes into that link -- either +-0.3 is a correct recovery).
gamma0    <- 0.3
prob0     <- exp(-abs(gamma0))
eff_ind   <- rbinom(n, 1, prob0)                ## 1 = efficient (u forced to 0)
y_zisf    <-         cons + beta1*x1    + beta2*x2   + v - u*(1-eff_ind)

## ZISF_Z -- targets zsfm(model_name="ZISF_Z", formula=y_zisf_z~x1+x2|z).
## ZISF_Z's prob link is the logistic prob = plogis(z'gamma); reusing the
## same z and 0.9/0.6 coefficients as the other z-parameterized columns
## above for consistency (even though here they parameterize a mixing
## probability, not a variance).
prob_z_true <- plogis(0.9 + 0.6*z)
eff_ind_z   <- rbinom(n, 1, prob_z_true)
y_zisf_z    <-         cons + beta1*x1    + beta2*x2   + v - u*(1-eff_ind_z)

## Output - sfm
y_pcs    <-         cons + beta1*x1    + beta2*x2   + v   - u
y_pcs_z  <-         cons + beta1*x1    + beta2*x2   + v   - uz
y_pcs_t  <-         cons + beta1*x1    + beta2*x2   + v_t - u_t
y_pcs_tn <-         cons + beta1*x1    + beta2*x2   + v   - u_tn
y_pcs_e  <-         cons + beta1*x1    + beta2*x2   + v   - u_e
y_pcs_ez <-         cons + beta1*x1    + beta2*x2   + v   - uz_e
y_pcs_c  <-         cons + beta1*x1    + beta2*x2   + v_c - u_c
y_pcs_u  <-         cons + beta1*x1    + beta2*x2   + v   - u_u

## Normal + Cauchy - half-normal
y_pcs_w  <-         cons + beta1*x1    + beta2*x2   + v +v_c   - u

data_trial <- as.data.frame(cbind(
  name,  cons,  x1,  x2,  u, uz,   v,  u_t,  v_t,  u_c,  v_c,  u_e , u_u,  u_tn,
  y_pcs,  y_pcs_t, y_pcs_e,   y_pcs_c,  y_pcs_u,   y_pcs_z,   y_pcs_w,   y_pcs_tn,   z,
  uz_e, y_pcs_ez,
  w_tt, y_ttne, w_tt_hn, y_tthn, zp, wz_hn, y_tthn_z,
  eff_ind, y_zisf, prob_z_true, eff_ind_z, y_zisf_z
))
colnames(data_trial) <- c(
  "name","cons","x1","x2","u","uz","v","u_t","v_t","u_c","v_c","u_e","u_u","u_tn",
  "y_pcs","y_pcs_t","y_pcs_e","y_pcs_c","y_pcs_u", "y_pcs_z", "y_pcs_w" ,"y_pcs_tn", "z",
  "uz_e", "y_pcs_ez",
  "w_tt", "y_ttne", "w_tt_hn", "y_tthn", "zp", "wz_hn", "y_tthn_z",
  "eff_ind", "y_zisf", "prob_z_true", "eff_ind_z", "y_zisf_z"
)

data_rand           <- data.frame(data_trial)
data_rand$con       <- 1
return(data_rand)
}
