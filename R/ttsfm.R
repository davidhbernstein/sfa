ttsfm <- function(formula, 
                model_name    = c("TTNE","TTHN"),
                data, 
                maxit.bobyqa  = 80000,
                maxit.psoptim = 1000,
                maxit.optim   = 1000, 
                REPORT        = 1, 
                trace         = 0,
                pgtol         = 0, 
                start_val     = FALSE,
                PSopt         = FALSE,
                optHessian    = TRUE,
                inefdec       = TRUE,
                upper         = NA,
                Method        = "L-BFGS-B",
                logit         = TRUE,
                verbose       = FALSE,
                rand.psoptim  = NULL){
  
call          <- match.call()
model_name    <- match.arg(model_name) 

DR1 <- data_proc(formula, data, model_name, individual = NULL, inefdec)

formula      <- DR1$formula
data_orig    <- DR1$data_orig
form_parts   <- DR1$form_parts
formula_x    <- DR1$formula_x
y_var        <- DR1$y_var
model_name   <- DR1$model_name
data_x       <- DR1$data_x
intercept    <- DR1$intercept
inefdec_n    <- DR1$inefdec_n
inefdec_TF   <- DR1$inefdec_TF
x_vars_vec   <- DR1$x_vars_vec
n_x_vars     <- DR1$n_x_vars
x_vars       <- DR1$x_vars
x_x_vec      <- DR1$x_x_vec
fancy_vars   <- DR1$fancy_vars
fancy_vars_z <- DR1$fancy_vars_z
n_z_vars     <- DR1$n_z_vars
N            <- DR1$N
data_z       <- DR1$data_z
if(length(unlist(form_parts))>3){  
  formula_z    <- DR1$formula_z
  intercept_z  <- DR1$intercept_z
  n_z_vars     <- DR1$n_z_vars
  z_vars       <- DR1$z_vars
  z_vars_vec   <- DR1$z_vars_vec
  z_z_vec      <- DR1$z_z_vec}
if(length(unlist(form_parts))>4){    ## might need to incorporate this above
  formula_zp    <- DR1$formula_zp
  intercept_zp  <- DR1$intercept_zp
  n_zp_vars     <- DR1$n_zp_vars
  zp_vars       <- DR1$zp_vars
  zp_vars_vec   <- DR1$zp_vars_vec
  zp_zp_vec     <- DR1$zp_zp_vec}


## Default starting values for variance equations
delta   <- rep(0.1, length(z_vars))
delta_p <- rep(0.1, length(zp_vars))
plm_lm      <- lm(formula_x ,data_orig)  
beta_0_st   <- if(isTRUE(intercept==0)) {NA} else{plm_lm$coefficients[c(1)]}
beta_hat    <- if(isTRUE(intercept==0)) {plm_lm$coefficients[x_vars_vec]} else{plm_lm$coefficients[x_vars_vec][-1]}
beta_0      <- beta_0_st  
sigma_v     <- .2

## Starting vector
if(isTRUE(is.numeric(start_val))){
  start_v <- start_val
} else {
  start_v <- if(is.na(beta_0_st)){
    unname(c(beta_hat, sigma_v, delta, delta_p))
  } else {
    unname(c(beta_0, beta_hat,sigma_v,delta, delta_p))
  }
}

## Output label matrix
out           <- matrix(0, nrow = 3, ncol = length(start_v))
rownames(out) <- c("par", "st_err", "t-val")
colnames(out) <- c(x_vars_vec,"sigv", z_vars, zp_vars)

DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num=NA, individual=NA, N, model_name, rand.gtre=NULL)

data         <- DR2$data
Y            <- DR2$Y

data_i_vars  <- DR2$data_i_vars
data_z_vars  <- as.matrix(data.frame(subset(data,select = z_vars)))
data_zp_vars <- as.matrix(data.frame(subset(data,select = zp_vars)))



if(model_name == "TTNE"){
fn = function(p){
    
  nr   <- n_x_vars ## number of regressors in regression
  nzu  <- n_z_vars ## number of determinants for u component
  nzw  <- n_zp_vars ## number of determinants for w component
  
  sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  sigu <- exp((data_z_vars%*%p[(nr+2):(nr+nzu+1)]))
  sigw <- exp((data_zp_vars%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))
  
  #if (sigv<= 1e-6){stop("Variance too small")}
  
  e <- Y - data_i_vars%*%p[1:nr]
  a <- (sigv^2)/(2*sigw^2) - e/sigw
  b <- e/sigv - sigv/sigw
  
  alpha <- e/sigu + (sigv^2)/(2*sigu^2)
  beta  <- -e/sigv - sigv/sigu
  
  denom <- sigu+sigw
  
  term1 <- exp(alpha)
  term2 <- exp(a)
  
  ##return will send the summation of the log of the 
  ##density of the composed error
  
  ll <- -log(denom)+log((pnorm(beta)*term1)+(pnorm(b)*term2))
  
  if(any(is.na(ll))){return(-.Machine$double.xmax)}
  if(is.null(ll)){return(-.Machine$double.xmax)}
  
  ll[ll==-Inf]         <-  -sqrt(.Machine$double.xmax/length(ll))
  ll[ll== Inf]         <-  -sqrt(.Machine$double.xmax/length(ll))
  ll[is.nan(ll)]       <-  -sqrt(.Machine$double.xmax/length(ll))
  
  return(sum(ll))
  
  # return(-sum(ll[is.finite(ll)]) )
}  
  
Start.Time <- start.time()

prep        <- list(n_x_vars, n_z_vars, n_zp_vars)
names(prep) <- c("n_x_vars", "n_z_vars", "n_zp_vars")

lower.BOB0 <- .generate_sfa_bounds(formula, prep)[-c(1:2)]
lower.BOB  <- append(lower.BOB0, .SFA_CONSTANTS$MIN_POSITIVE, after = n_x_vars)


## ---- Stage 1: bobyqa
Opt.Bobyqa <- opt.bobyqa(
  fn = fn,
  start_v = start_v,
  lower.bobyqa = lower.BOB,
  maxit.bobyqa = maxit.bobyqa,
  bob.TF = TRUE,
  verbose = verbose
)

start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1

## ---- Stage 2: psoptim
differ <- 10
lower1_0 <- .generate_sfa_bounds(formula, prep , inf_sub = min(start_v[-c(n_x_vars+1)])-differ )[-c(1:2)]
lower1   <- append(lower1_0, .SFA_CONSTANTS$MIN_POSITIVE, after = n_x_vars)   

Opt.Psoptim <- opt.psoptim(
  fn = fn,
  start_v,
  lower.psoptim = lower1,
  upper.psoptim = c(start_v + differ),
  rand.psoptim = rand.psoptim,
  maxit.psoptim = maxit.psoptim,
  psopt.TF = PSopt,
  verbose = verbose
)

start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

## ---- Stage 3: optim
differ <- 1
lower1_0 <- .generate_sfa_bounds(formula, prep , inf_sub = min(start_v[-c(n_x_vars+1)])-differ )[-c(1:2)]
lower1   <- append(lower1_0, .SFA_CONSTANTS$MIN_POSITIVE, after = n_x_vars)   

Opt.Optim <- opt.optim(
  fn = fn,
  start_v = start_v,
  lower.optim = lower1,
  upper.optim = c(start_v + differ),
  maxit.optim = maxit.optim,
  opt.TF = optHessian,
  method = Method,
  optHessian = optHessian,
  verbose = verbose
)

start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)

## Preserve current fallback logic
if(optHessian == FALSE && PSopt == FALSE){
  opt <- bob1
}

if(optHessian == FALSE && PSopt == TRUE){
  opt <- opt00
}


## now for st errs 
if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

if(optHessian==TRUE){ st_err  <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0 ) ){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(diag(solve(opt$hessian))))}}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val

## metrics 
metrics.ne <- function(p, e=NULL, y=Y, xx=data_i_vars, zu=data_z_vars, zw=data_zp_vars, alphahat=NULL){
  
  if(is.null(e)){
    
    nr  <- ncol(xx)  ##Calculate number of regressors in regression
    nzu <- ncol(zu)  ##Calculate number of determinants for u component
    nzw <- ncol(zw) ##Calculate number of determinants for w component
    
    ep.hat  <- y-xx%*%p[1:nr]
    
    if(!is.null(alphahat)){
      
      ep.hat  <- y-xx%*%p[1:nr]-alpha.hat
      
    }
    
    
    sig.v <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
    sig.u <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
    sig.w <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))
    
  }else{
    
    ep.hat  <- e
    
    sig.v  <- exp(p[1]) 	
    sig.u  <- exp(p[2])
    sig.w <- exp(p[3])
    
    ## Have to correct for the shift in our residuals to begin with
    ep.hat <- e-sig.u+sig.w
    
  }
  ## Use 8.23 and 8.26 to construct metrics
  ## Setup necessary parameters needed.
  lambda <- 1/sig.w+1/sig.u
  a1     <- sig.v^2/(2*sig.u^2)+ep.hat/sig.u
  b1     <- -(ep.hat/sig.v+sig.v/sig.u)
  a2     <- sig.v^2/(2*sig.w^2)-ep.hat/sig.w
  b2     <-  ep.hat/sig.v-sig.v/sig.w
  chi1   <- pnorm(b2)+exp(a1-a2)*pnorm(b1)
  chi2   <- exp(a2-a1)*chi1
  
  Eew.cond  <- (lambda/(chi2*(lambda-1)))*(pnorm(b1)+
                                             exp(0.5*((b2+sig.v)^2-b1^2))*pnorm(b2+sig.v))
  
  Eemw.cond <- (lambda/(chi2*(1+lambda)))*(pnorm(b1)+
                                             exp(a2-a1-b2*sig.v+0.5*sig.v^2)*pnorm(b2-sig.v))
  
  Eeu.cond  <- (lambda/(chi1*(lambda-1)))*(pnorm(b2)+
                                             exp(0.5*((b1+sig.v)^2-b2^2))*pnorm(b1+sig.v))
  
  Eemu.cond <- (lambda/(chi1*(1+lambda)))*(pnorm(b2)+
                                             exp(a1-a2-b1*sig.v+0.5*sig.v)*pnorm(b1-sig.v))
  
  Eewmu.cond <- (exp((1+sig.u)*(a1+sig.v^2/2/sig.u))*pnorm(b1-sig.v)+
                   exp((1-sig.w)*(a2-sig.v^2/2/sig.w))*
                   pnorm(b2+sig.v))/(exp(a1)*pnorm(b1)+
                                       exp(a2)*pnorm(b2))
  
  Eeumw.cond <- (exp((1-sig.u)*(a1-sig.v^2/2/sig.u))*pnorm(b1+sig.v)+
                   exp((1+sig.w)*(a2+sig.v^2/2/sig.w))*
                   pnorm(b2-sig.v))/(exp(a1)*pnorm(b1)+
                                       exp(a2)*pnorm(b2))
  
  ## Now calculate the M1 and M2 metrics (these are information deficiency 
  ## relative to the actual price) and M5 and M6 metrics (these are information 
  ## deficiency relative to balanced price).
  M1.ne <- 1-Eemw.cond
  M2.ne <- Eeu.cond-1
  
  M5.ne <- Eew.cond-1
  M6.ne <- 1-Eemu.cond
  
  ## Alecos wanted M7 as well
  M7.ne   <- Eewmu.cond-1
  M10.ne <- 1-Eeumw.cond
  
  return(list(M1=M1.ne, M2=M2.ne, M5=M5.ne,M6=M6.ne, M7=M7.ne, M10=M10.ne, Eew.cond=Eew.cond,Eemw.cond=Eemw.cond,
              Eeu.cond=Eeu.cond, Eemu.cond=Eemu.cond, Eewmu.cond=Eewmu.cond, Eeumw.cond=Eeumw.cond))
  
}

metric.ne.res <- metrics.ne(p=opt$par)


results <- list(t(out),c(opt),End.Time,start_v,model_name,formula,out["par",], out["st_err",], out["t-val",],call)
class(results)  <- "sfareg"
names(results)  <- c("out","opt","total_time","start_v","model_name","formula","coefficients", "std.errors", "t.values","call")
return(results)}

else{return(c("This is not a valid command"))}}



