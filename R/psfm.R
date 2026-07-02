psfm <- function(formula, 
                 model_name   = c("TRE_Z","GTRE_Z","TRE","GTRE","TFE","FD","GTRE_SEQ1","GTRE_SEQ2"), 
                 data, 
                 maxit.bobyqa = 100,
                 maxit.psoptim= 10,
                 maxit.optim  = 10,
                 REPORT       = 1,
                 trace        = 3, 
                 pgtol        = 0, 
                 individual,
                 halton_num   = NULL,
                 start_val    = FALSE,
                 gamma        = FALSE,
                 PSopt        = FALSE,
                 optHessian   = TRUE,
                 inefdec      = TRUE, 
                 Method       = "L-BFGS-B",
                 verbose      = FALSE,
                 rand.gtre    = NULL,
                 rand.psoptim = NULL){
    
call          <- match.call()
model_name    <- match.arg(model_name)
  
DR1 <- data_proc(formula, data, model_name, individual, inefdec)

if(1==1){
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

Start_Panel  <- start_panel(formula_x, data, model_name, start_val, intercept, x_vars_vec)

alpha_hat    <- Start_Panel$alpha_hat
beta_0       <- Start_Panel$beta_0
beta_0_st    <- Start_Panel$beta_0_st
beta_hat     <- Start_Panel$beta_hat
beta_se      <- Start_Panel$beta_se
epsilon_hat  <- Start_Panel$epsilon_hat
exp_eta      <- Start_Panel$exp_eta
exp_u        <- Start_Panel$exp_u
lambda       <- Start_Panel$lambda
out          <- Start_Panel$out
plm_gtre     <- Start_Panel$plm_gtre
sfa_alp      <- Start_Panel$sfa_alp
sfa_eps      <- Start_Panel$sfa_eps
sigma        <- Start_Panel$sigma
sigma_h      <- Start_Panel$sigma_h
sigma_r      <- Start_Panel$sigma_r
sigma_u      <- Start_Panel$sigma_u
sigma_v      <- Start_Panel$sigma_v
start_v      <- Start_Panel$start_v
plm_tfe      <- Start_Panel$plm_tfe
plm_fd       <- Start_Panel$plm_fd

DR2 <- data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, 
                  x_vars_vec, halton_num, individual, N, model_name, rand.gtre)

data         <- DR2$data
Y            <- DR2$Y
data_i_vars  <- DR2$data_i_vars}

if(model_name %in% c("GTRE","TRE")){
  
R      <- DR2$R
R_H    <- DR2$R_H
indiv  <- DR2$indiv
t      <- DR2$t
data_i <- DR2$data_i
eps    <- DR2$eps 
R_h1   <- DR2$R_h1
R_h2   <- DR2$R_h2
data_x <- DR2$data_x
  
fn_1 = function(x){
   if(model_name == "GTRE"){x_x_vec       <- x[5:as.numeric(n_x_vars + 4)]}
   if(model_name == "TRE"){ x_x_vec       <- x[4:as.numeric(n_x_vars + 3)]}
 fn1 = function(ii){
   if(model_name == "GTRE"){eps_neg    <- eps 
   eps_neg[[ii]]     <- Y[[ii]]  + x[3]*R_h1[[ii]] + x[4]*R_h2[[ii]] * inefdec_n}
   
   if(model_name == "GTRE"){eps[[ii]]         <- Y[[ii]]  - x[3]*R_h1[[ii]] + x[4]*R_h2[[ii]] * inefdec_n}
   if(model_name == "TRE"){ eps[[ii]]         <- Y[[ii]]  - x[3]*R_h1[[ii]]}
    
    for (qq in 1:n_x_vars){
    eps[[ii]] <- eps[[ii]] - x_x_vec[qq]*matrix(rep(data_i_vars[[ii]][,qq],R),t[[ii]],R)}
   
   if(model_name == "GTRE"){ 
     for (qq in 1:n_x_vars){
       eps_neg[[ii]] <- eps_neg[[ii]] - x_x_vec[qq]*matrix(rep(data_i_vars[[ii]][,qq],R),t[[ii]],R)} 
     
     eps_neg[[ii]] <- inefdec_n*eps_neg[[ii]]
     z1_neg        <-  eps_neg[[ii]]/x[2]
     z2_neg        <- -eps_neg[[ii]]*x[1]/x[2]} 
   
    eps[[ii]] <- inefdec_n*eps[[ii]]
    z0          <-  2/x[2]
    z1          <-  eps[[ii]]/x[2]
    z2          <- -eps[[ii]]*x[1]/x[2]
    z1[z1>  .SFA_CONSTANTS$CLIP_Z1_UPPER] <- .SFA_CONSTANTS$CLIP_Z1_UPPER 
    z1[z1<  .SFA_CONSTANTS$CLIP_Z1_LOWER] <- .SFA_CONSTANTS$CLIP_Z1_LOWER
    z2[z2>  .SFA_CONSTANTS$CLIP_Z2_UPPER] <- .SFA_CONSTANTS$CLIP_Z2_UPPER
    z2[z2<  .SFA_CONSTANTS$CLIP_Z2_LOWER] <- .SFA_CONSTANTS$CLIP_Z2_LOWER
            
            prod_vec_n0  <- log(max(   mean(   colProds(    
                         z0* dnorm(z1)* pmax(pnorm(z2), eps[[ii]]*0+.Machine$double.xmin)    
                       )), .Machine$double.xmin))
            if(model_name == "TRE"){ prod_vec_n <- prod_vec_n0}
            if(model_name == "GTRE"){ 
            prod_vec_n1  <- log(max(   mean(   colProds(    
              z0* dnorm(z1_neg)* pmax(pnorm(z2_neg), eps_neg[[ii]]*0+.Machine$double.xmin)    
            )), .Machine$double.xmin))
            
            prod_vec_n <- 0.5 * (prod_vec_n0 + prod_vec_n1) }
            
            return(-prod_vec_n)}
          
fn1_apply                         <- unlist(lapply(1:N, fn1))
fn1_apply[is.nan(fn1_apply)]      <- sqrt(.SFA_CONSTANTS$MAX_VALUE/length(x))
fn1_apply[is.infinite(fn1_apply)] <- sqrt(.SFA_CONSTANTS$MAX_VALUE/length(x))
          
return( sum( fn1_apply[is.finite(fn1_apply)] ) ) }

Start.Time  <- start.time()      
Lower.Start <- lower.start(start_v, model_name, differ=3)
Opt.Bobyqa  <- opt.bobyqa(fn=fn_1, start_v=start_v, lower.bobyqa=Lower.Start$lower1, maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose)
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 
      
Lower.Start <- lower.start(start_v, model_name, differ=2)
Opt.Psoptim <- opt.psoptim(fn=fn_1, start_v=start_v, lower.psoptim=Lower.Start$lower1,
      rand.psoptim=rand.psoptim,upper.psoptim=Lower.Start$upper1, maxit.psoptim, psopt.TF=PSopt,
      verbose = verbose)
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

Lower.Start <-lower.start(start_v, model_name, differ=0.5)
Opt.Optim <- opt.optim(fn = fn_1, start_v = start_v, lower.optim =Lower.Start$lower1 ,upper.optim=Lower.Start$upper1, 
          maxit.optim=maxit.optim, opt.TF = optHessian, method=Method, optHessian= optHessian,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)         

if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))} 
out            <- matrix(0,nrow = 3,ncol = length(start_v))
rownames(out)  <- c("par","st_err","t-val")
colnames(out)  <- if(model_name=="GTRE"){c("lambda","sigma","sigr","sigh",colnames(data_x))} else{c("lambda","sigma","sigr",colnames(data_x))} 

st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE ){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)),0)))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
      
## TE Measurements : GTRE       
if(model_name == "GTRE"){ 
        
        beta  <- opt$par[-c(1:4)]
        lamb  <- opt$par[1]
        sig   <- opt$par[2]
        sig_u <- (lamb*sig) / sqrt(1+lamb^2)
        sig_v <- sig_u/lamb
        sig_r <- opt$par[3]
        sig_h <- opt$par[4]
        
        e_i  <-   as.list(rep(0,N))
        A    <-   as.list(rep(0,N))
        SIG  <-   as.list(rep(0,N))
        VEE  <-   as.list(rep(0,N))
        LAM  <-   as.list(rep(0,N))
        ARR  <-   as.list(rep(0,N))
        n    <-   sum(t)
        U    <-   rep(0,n)
        
        for(i in 1:N){
          e_i[[i]]     <- pmin(inefdec_n*(Y[[i]][,1] - rowSums(t(t(data_i_vars[[i]]) * beta))),Y[[i]][,1]*0)
          A[[i]]       <- -cbind(rep(1,t[i]),diag(t[i]))
          SIG[[i]]     <- sig_v^2*diag(t[i]) + sig_r^2*rep(1,t[i])%*%t(rep(1,t[i]))
          VEE[[i]]     <- rbind( c(sig_h^2,rep(0,t[i]))   ,  cbind( rep(0,t[i]) , sig_u^2*diag(t[i])) )
          LAM[[i]]     <- solve( solve(VEE[[i]])  +  t(A[[i]]) %*% solve(SIG[[i]]) %*% A[[i]]    )
          ARR[[i]]     <- LAM[[i]] %*% t(A[[i]]) %*% solve(SIG[[i]])
        }
        
        res_d_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx = rep(Inf, t[i]+1), 
                                         mean = as.numeric(ARR[[i]] %*% e_i[[i]] ),
                                         sigma=LAM[[i]])[1]}
        
        res_n_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx=rep(Inf, t[i]+1), 
                                         mean =  as.numeric(ARR[[i]] %*% e_i[[i]] + LAM[[i]]%*% c(-1,rep(0,t[i])) ), 
                                         sigma=LAM[[i]])[1]}
        
        res_d <- lapply(seq(1,N,1), res_d_fn)
        res_n <- lapply(seq(1,N,1), res_n_fn)
        
        H_fn    <- function(i){(max(res_n[[i]], .SFA_CONSTANTS$MIN_POSITIVE)/max(res_d[[i]],.SFA_CONSTANTS$MIN_POSITIVE) )*
            exp(t(c(-1,rep(0,t[i])))%*%ARR[[i]]%*%e_i[[i]] + 0.5*t(c(-1,rep(0,t[i])))%*%LAM[[i]]%*%c(-1,rep(0,t[i])) )}
        
        H    <- unlist(lapply(seq(1,N,1), H_fn))
        H    <- pmin(H,rep(1,length(H)))
        
        new_t_exp <- as.list(rep(0,n))
        
        for(i in 1:N){
          for (j in 1:t[i]) {
            h <- cumsum(t)[i]-t[i] + j
            new_t      <- rep(0,t[i] +1)
            new_t[j+1] <- -1 
            new_t_exp[[h]] <- new_t
            
          }}
        
        
t_cum      <-   c(cumsum(t))
t_exp      <-   e_i_exp    <-   A_exp  <- SIG_exp <- VEE_exp <- LAM_exp <- ARR_exp <- res_d_exp  <-   as.list(rep(0,n))
        
for(m in 1:N){
B  <- t_cum[m]
A  <- B +1 - t[m] 
t_exp[A:B]     <- rep(t[m],  t[m])
e_i_exp[A:B]   <- rep(e_i[m],t[m])
A_exp[A:B]     <- rep(A[m],  t[m])
SIG_exp[A:B]   <- rep(SIG[m],t[m])
VEE_exp[A:B]   <- rep(VEE[m],t[m])
LAM_exp[A:B]   <- rep(LAM[m],t[m])
ARR_exp[A:B]   <- rep(ARR[m],t[m])
res_d_exp[A:B] <- rep(res_d[m],t[m])
}
        
        res_n_t_fn <- function(i){ptmvnorm(  lowerx=rep(0, t_exp[[i]]+1), upperx=rep(Inf, t_exp[[i]]+1),
                                             mean= as.numeric(ARR_exp[[i]] %*% e_i_exp[[i]] + LAM_exp[[i]]%*% new_t_exp[[i]] ), 
                                             sigma=LAM_exp[[i]])[1]}
        
        res_n_t <- lapply(seq(1,n,1), res_n_t_fn)
        
        
        U_fn    <- function(i){(max(res_n_t[[i]], .SFA_CONSTANTS$MIN_POSITIVE) /  max(res_d_exp[[i]], .SFA_CONSTANTS$MIN_POSITIVE))*
            exp(t(new_t_exp[[i]])%*%ARR_exp[[i]]%*%e_i_exp[[i]] + 0.5*t(new_t_exp[[i]])%*%LAM_exp[[i]]%*%new_t_exp[[i]] )}
        
        U     <- unlist(lapply(seq(1,n,1), U_fn))
        U     <- pmin(U,rep(1,length(U)))}
## TE Measurements : TRE       
if(model_name == "TRE"){
        beta  <- opt$par[-c(1:3)]
        lamb  <- opt$par[1]
        sig   <- opt$par[2]
        sig_u <- (lamb*sig) / sqrt(1+lamb^2)
        sig_v <- sig_u/lamb
        
        Y_mean           <- rep(0,N)
        X_mean           <- matrix(0,N,ncol = length(x_vars_vec))
        colnames(X_mean) <- x_vars_vec
        
        for (ii in 1:N) {
          data_i[[ii]]  <- data[which(data[,c(individual)]==indiv[ii]),]
          Y_mean[ii]    <- mean(as.numeric(data_i[[ii]][,y_var]))
          X_mean[ii,]   <- colMeans(data.frame(data_i[[ii]][,c(x_vars_vec)] ))
        }
        
        r_hat_m      <- Y_mean - rowSums(t(beta*t(X_mean))) + inefdec_n*sqrt(2/pi)*sig_u
        r_hat_m_exp  <- rep(0,sum(t))
        t_cum        <- c(cumsum(t))
        
        for(m in 1: length(t)){
          B  <- t_cum[m]
          A  <- B +1 - t[m] 
          r_hat_m_exp[A:B] <- rep(r_hat_m[m],t[m])}
        
eps_hat    <- pmin(inefdec_n*(data[,y_var] - rowSums(t(t(data[,c(x_vars_vec)])*beta)) -r_hat_m_exp ),data[,y_var]*0)
sig_star   <- sig_u*sig_v/sig
inner      <- (lamb*eps_hat)/sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/(1-pnorm(inner)))*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
U          <- exp_u_hat}
      
#cor(U,exp(-p_data_trial$u))
#cor(H,exp(-unique(p_data_trial$h) ))

if(model_name == "GTRE"){      
results <- list(t(out),c(opt),data,End.Time,start_v,model_name,formula, U, H, out["par",], out["st_err",], out["t-val",], call)
class(results)  <- "sfareg"
names(results)  <- c("out","opt","data","total_time","start_v","model_name","formula","U","H", "coefficients", "std.errors", "t.values","call")}
if(model_name == "TRE"){      
results <- list(t(out),c(opt),data,End.Time,start_v,model_name,formula, U, out["par",], out["st_err",], out["t-val",], call)
class(results)  <- "sfareg"
names(results)  <- c("out","opt","data","total_time","start_v","model_name","formula","U", "coefficients", "std.errors", "t.values", "call")}

return(results)}   
if(model_name ==     "GTRE_Z")     { 
## Default starting values for variance equations
    delta   <- rep(0.1, length(z_vars))
    delta_p <- rep(0.1, length(zp_vars))
    
    ## Starting vector
    if(isTRUE(is.numeric(start_val))){
      start_v <- start_val
    } else {
      start_v <- if(is.na(beta_0_st)){
        unname(c(sigma_v, sigma_r, beta_hat, delta, delta_p))
      } else {
        unname(c(sigma_v, sigma_r, beta_0, beta_hat, delta, delta_p))
      }
    }
    
    ## Output label matrix
    out           <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- c("sigv", "sigr", colnames(data[, x_vars_vec, drop = FALSE]), z_vars, zp_vars)
    
    ## Number of simulation draws
    R <- if(isTRUE(is.numeric(halton_num))) halton_num else ceiling(sqrt(nrow(data))) + 100
    
    ## Halton draws
    R_H <- randtoolbox::halton(
      R + .SFA_CONSTANTS$HALTON_DISCARD,
      2,
      start = 1,
      normal = FALSE)[-c(1:.SFA_CONSTANTS$HALTON_DISCARD), c(1:2)]
    
    ## First column: standard normal via qnorm
    ## Second column: half-normal via inverse error function
    R_H <- cbind(qnorm(R_H[,1]), sqrt(2) * pracma::erfinv(R_H[,2]))
    
    if(!is.null(rand.gtre)){
      set.seed(rand.gtre)}
    
    ## Optional decorrelation step retained from current code
    mat <- matrix(0, nrow = R, ncol = 9999)
    for(v in 1:9999){
      mat[,v] <- sample(R_H[,1])
    }
    
    cor_vec <- matrix(0, 9999, 1)
    for(v in 1:9999){
      cor_vec[v] <- abs(cor(mat[,v], R_H[,2]))
    }
    
    R_H <- cbind(mat[, which(cor_vec == min(cor_vec))[1]], R_H[,2])
    rm(mat, cor_vec, v)
    
    ## Build firm-level lists
    indiv <- noquote(as.vector(unique(data[, individual])))
    t_vec <- rep(0, N)
    
    data_i       <- vector("list", N)
    Y            <- vector("list", N)
    data_i_vars  <- vector("list", N)
    data_z_vars  <- vector("list", N)
    data_zp_vars <- vector("list", N)
    R_h1         <- vector("list", N)
    R_h2         <- vector("list", N)
    
    for(ii in seq_len(N)){
      data_i[[ii]]       <- data[which(data[, individual] == indiv[ii]), , drop = FALSE]
      t_vec[ii]          <- nrow(data_i[[ii]])
      R_h1[[ii]]         <- t(matrix(rep(R_H[,1], t_vec[[ii]]), R, t_vec[[ii]]))
      R_h2[[ii]]         <- abs(t(matrix(rep(R_H[,2], t_vec[[ii]]), R, t_vec[[ii]])))
      Y[[ii]]            <- matrix(rep(data_i[[ii]][, y_var], R), t_vec[[ii]], R)
      data_i_vars[[ii]]  <- data.frame(data_i[[ii]][, x_vars_vec, drop = FALSE])
      data_z_vars[[ii]]  <- data.frame(data_i[[ii]][, z_vars, drop = FALSE])
      data_zp_vars[[ii]] <- data.frame(data_i[[ii]][, zp_vars, drop = FALSE])
    }
    
    prep <-  list(
      data         = data,                  n_z_vars     = length(z_vars),   t            = t_vec,
      individual   = individual,            n_zp_vars    = length(zp_vars),  data_i       = data_i,
      y_var        = y_var,                 N            = N,                Y            = Y,
      x_vars_vec   = x_vars_vec,            R            = R,                data_i_vars  = data_i_vars,
      z_vars       = z_vars,                start_v      = start_v,          data_z_vars  = data_z_vars,
      zp_vars      = zp_vars,               out_template = out,              data_zp_vars = data_zp_vars,
      n_x_vars     = n_x_vars,              indiv        = indiv,            R_h1         = R_h1,
      R_h2         = R_h2,
      formula      = formula)
  
    t            <-  t_vec  
    n_z_vars     <-  length(z_vars)
    n_zp_vars    <-  length(zp_vars)
    out_template <-  out
  
    ## GTRE simulated negative log-likelihood
    fn <-  function(x){
      
      beta_start   <- 3
      beta_end     <- beta_start + n_x_vars - 1
      
      delta_start  <- beta_end + 1
      delta_end    <- delta_start + n_z_vars - 1
      
      deltap_start <- delta_end + 1
      deltap_end   <- deltap_start + n_zp_vars - 1
      
      beta         <- x[beta_start:beta_end]
      
      delta        <- x[delta_start:delta_end]
      delta_p      <- x[deltap_start:deltap_end]
      
      sig_v        <- x[1]
      sig_r        <- x[2]
      
      ## Enforce positivity on sigma_v and sigma_r
      sig_v <- max(sig_v, .SFA_CONSTANTS$MIN_POSITIVE)
      sig_r <- max(sig_r, .SFA_CONSTANTS$MIN_POSITIVE)
      
      ## Bound the variance linear predictors before exponentiating.
      eta_bound <- 40
      
      ## One firm-level log-likelihood contribution
      ll_i <- function(ii){
        
        Ti <- t[ii]
        
        ## Persistent inefficiency scale for firm i
        ## IMPORTANT:
        ## This is computed inside the firm-specific contribution,
        ## not outside, so it is truly unit-specific.
        sigma_h_fun <- if(n_zp_vars > 1){
          eta_h_i <- as.numeric(as.matrix(data_zp_vars[[ii]]) %*% delta_p)
          if(any(!is.finite(eta_h_i))){
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)}
          eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
          mean(sqrt(exp(eta_h_i)))
        } else {
          eta_h_i <- x[deltap_end]
          if(!is.finite(eta_h_i)){return((.SFA_CONSTANTS$MAX_VALUE)^0.1)}
          eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
          sqrt(exp(eta_h_i))}
        
        
        ## Construct epsilon_it draw-by-draw
        eps_ii <- Y[[ii]] - sig_r * R_h1[[ii]] + sigma_h_fun * R_h2[[ii]] * inefdec_n
        
        ## Remove frontier mean
        for(qq in seq_len(n_x_vars)){
          eps_ii <- eps_ii - beta[qq] *
            matrix(rep(data_i_vars[[ii]][, qq], R), Ti, R)
        }
        
        eps_ii <- inefdec_n * eps_ii
        
        ## Transient inefficiency scale
        sigma_u_fun <- if(n_z_vars > 1){
          eta_u_i <- as.numeric(as.matrix(data_z_vars[[ii]]) %*% delta)
          if(any(!is.finite(eta_u_i))){
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
          }
          eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
          sqrt(exp(eta_u_i))
        } else {
          eta_u_i <- x[delta_end]
          if(!is.finite(eta_u_i)){return((.SFA_CONSTANTS$MAX_VALUE)^0.1)}
          eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
          sqrt(exp(eta_u_i))
        }
        
        sigma_fun <- sqrt(sig_v^2 + sigma_u_fun^2)
        lamb_fun  <- sigma_u_fun / sig_v
        
        ## Expand firm-level sigma objects across draws
        sigma_mat <- matrix(rep(sigma_fun, R), Ti, R)
        lamb_mat  <- matrix(rep(lamb_fun,  R), Ti, R)
        
        ## Simulated likelihood contribution:
        ##   average over R draws of the product across t
        sim_terms <- (2 / sigma_mat) *
          dnorm(eps_ii / sigma_mat) *
          pmax(
            pnorm(-eps_ii * lamb_mat / sigma_mat),
            eps_ii * 0 + .SFA_CONSTANTS$MIN_POSITIVE)
        prod_vec_n <- log(mean(colProds(sim_terms)))
        -prod_vec_n
      }
      
      ll_vec <- unlist(lapply(seq_len(N), ll_i))
      
      ll_vec[which(ll_vec == Inf)]  <-  (.SFA_CONSTANTS$MAX_VALUE)^0.1
      ll_vec[which(ll_vec == -Inf)] <- -(.SFA_CONSTANTS$MAX_VALUE)^0.1
      
      sum(ll_vec[is.finite(ll_vec)])
    }
    
  ## Staged optimization for GTRE
    Start.Time <- start.time()
    
    lower.BOB <- .generate_sfa_bounds(formula, prep)  # default to -Inf 
    
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
    # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
    lower1 <- .generate_sfa_bounds(formula, prep , inf_sub = min(start_v[-c(1:2)])-differ )   
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
    # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
    lower1 <- .generate_sfa_bounds(formula, prep , inf_sub = min(start_v[-c(1:2)])-differ )   
    
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
      opt <- bob1}
    
    if(optHessian == FALSE && PSopt == TRUE){
      opt <- opt00}
    
  
  
  ## ------------------------------------------------------------
  ## Post-estimation GTRE technical efficiency recovery
  ##
  ## This uses the refactored TE code we developed earlier.
  ## ------------------------------------------------------------
  .gtre_te <- function(opt, prep, inefdec_n){
    
    ## Basic dimensions and indexing
    n       <- sum(t)
    id_obs  <- rep(seq_len(N), t)
    t_cum   <- cumsum(t)
    t_start <- c(1, head(t_cum, -1) + 1)
    
    ## Build variance-design matrices
    Z_mat  <- .make_var_design(data, z_vars,  rows = NULL,    int_name = "int_u")
    Zp_mat <- .make_var_design(data, zp_vars, rows = t_start, int_name = "int_h")
    
    # n_z_eff  <- ncol(Z_mat)
    # n_zp_eff <- ncol(Zp_mat)
    
    ## Parameter indexing
    beta_start   <- 3
    beta_end     <- beta_start + n_x_vars - 1
    
    delta_start  <- beta_end + 1
    delta_end    <- delta_start + n_z_vars - 1
    
    deltap_start <- delta_end + 1
    deltap_end   <- deltap_start + n_zp_vars - 1
    
    beta    <- opt$par[beta_start:beta_end]
    delta   <- opt$par[delta_start:delta_end]
    delta_p <- opt$par[deltap_start:deltap_end]
    
    sig_v <- max(opt$par[1], 1e-8)
    sig_r <- max(opt$par[2], 1e-8)
    
    min_sd    <- 1e-8
    eta_bound <- 40
    
    ## Compute and check the variance linear predictors before exp().
    ## The earlier pmax() guard only protected against near-zero standard
    ## deviations.  It did not protect against overflow from exp(Z %*% delta),
    ## which can produce Inf values and cause the posterior covariance system
    ## to fail later in .safe_inverse().
    eta_u <- as.numeric(Z_mat  %*% delta)
    eta_h <- as.numeric(Zp_mat %*% delta_p)
    
    if(any(!is.finite(eta_u)) || any(!is.finite(eta_h))){
      stop(
        "Non-finite variance linear predictor in .gtre_te(). Check Z*delta and Zp*delta_p.",
        call. = FALSE
      )
    }
    
    ## Bound before exponentiating to prevent numerical overflow.
    eta_u <- pmin(pmax(eta_u, -eta_bound), eta_bound)
    eta_h <- pmin(pmax(eta_h, -eta_bound), eta_bound)
    
    sig_u_all <- pmax(sqrt(exp(eta_u)), min_sd)
    sig_h_all <- pmax(sqrt(exp(eta_h)), min_sd)
    
    if(any(!is.finite(sig_u_all)) || any(!is.finite(sig_h_all))){
      stop(
        "Non-finite sigma_u or sigma_h values in .gtre_te(). Check bounded variance predictors.",
        call. = FALSE
      )
    }
    
    sig_u_split <- split(sig_u_all, id_obs)
    
    ## Residual vectors epsilon_i
    e_i <- Map(
      f = function(Yi, Xi){
        pmin(
          inefdec_n * (Yi[,1] - rowSums(t(t(Xi) * beta))),
          Yi[,1] * 0
        )
      },
      Yi = Y,
      Xi = data_i_vars
    )
    
    ## Firm-level matrices
    A_i <- lapply(t, function(Ti) -cbind(rep(1, Ti), diag(Ti)))
    
    SIG <- lapply(
      t,
      function(Ti) sig_v^2 * diag(Ti) + sig_r^2 * tcrossprod(rep(1, Ti))
    )
    
    VEE <- Map(
      f = function(sig_h_i, sig_u_i){
        Ti <- length(sig_u_i)
        rbind(
          c(sig_h_i^2, rep(0, Ti)),
          cbind(rep(0, Ti), diag(sig_u_i^2, nrow = Ti, ncol = Ti))
        )
      },
      sig_h_i = sig_h_all,
      sig_u_i = sig_u_split
    )
    
    ## Build posterior covariance objects one firm at a time so that
    ## any numerical failure reports the offending firm and panel length.
    ## This makes failures much easier to debug than a generic Map/mapply
    ## traceback.
    post_obj <- lapply(seq_len(N), function(ii){
      if(any(!is.finite(VEE[[ii]])) || any(!is.finite(A_i[[ii]])) || any(!is.finite(SIG[[ii]]))){
        stop(
          sprintf(
            "Non-finite input matrix in .gtre_te() for firm %s (Ti = %s).",
            ii, t[ii]
          ),
          call. = FALSE
        )
      }
      
      tryCatch(
        .safe_linear_combo(
          VEE_i      = VEE[[ii]],
          A_i        = A_i[[ii]],
          SIG_i      = SIG[[ii]],
          base_ridge = 1e-10,
          ridge_mult = 10,
          max_tries  = 8
        ),
        error = function(e){
          stop(
            sprintf(
              "Failed to compute GTRE posterior matrices in .gtre_te() for firm %s (Ti = %s): %s",
              ii, t[ii], conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
    })
    
    LAM <- lapply(post_obj, `[[`, "LAM")
    ARR <- lapply(post_obj, `[[`, "ARR")
    
    ridge_report <- data.frame(
      firm        = seq_len(N),
      ridge_VEE   = sapply(post_obj, `[[`, "ridge_VEE"),
      ridge_SIG   = sapply(post_obj, `[[`, "ridge_SIG"),
      ridge_K     = sapply(post_obj, `[[`, "ridge_K"),
      method_VEE  = sapply(post_obj, `[[`, "invVEE_method"),
      method_SIG  = sapply(post_obj, `[[`, "invSIG_method"),
      method_K    = sapply(post_obj, `[[`, "invK_method")
    )
    
    ## Persistent TE
    res_d <- mapply(
      FUN = function(Ti, ARR_i, e_i, LAM_i){
        ptmvnorm(
          lowerx = rep(0, Ti + 1),
          upperx = rep(Inf, Ti + 1),
          mean   = as.numeric(ARR_i %*% e_i),
          sigma  = LAM_i
        )[1]
      },
      Ti       = t,
      ARR_i    = ARR,
      e_i      = e_i,
      LAM_i    = LAM,
      SIMPLIFY = FALSE
    )
    
    res_n <- mapply(
      FUN = function(Ti, ARR_i, e_i, LAM_i){
        shift_vec <- c(-1, rep(0, Ti))
        ptmvnorm(
          lowerx = rep(0, Ti + 1),
          upperx = rep(Inf, Ti + 1),
          mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
          sigma  = LAM_i
        )[1]
      },
      Ti       = t,
      ARR_i    = ARR,
      e_i      = e_i,
      LAM_i    = LAM,
      SIMPLIFY = FALSE
    )
    
    H <- mapply(
      FUN = function(Ti, ARR_i, e_i, LAM_i, rd, rn){
        shift_vec <- c(-1, rep(0, Ti))
        (max(rn, .SFA_CONSTANTS$MIN_POSITIVE) /
            max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
          exp(
            t(shift_vec) %*% ARR_i %*% e_i +
              0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
          )
      },
      Ti    = t,
      ARR_i = ARR,
      e_i   = e_i,
      LAM_i = LAM,
      rd    = res_d,
      rn    = res_n
    )
    
    H <- pmin(H, 1)
    
    ## Transient TE
    U_list <- mapply(
      FUN = function(Ti, ARR_i, e_i, LAM_i, rd){
        
        sapply(seq_len(Ti), function(j){
          
          shift_vec <- rep(0, Ti + 1)
          shift_vec[j + 1] <- -1
          
          rn_t <- ptmvnorm(
            lowerx = rep(0, Ti + 1),
            upperx = rep(Inf, Ti + 1),
            mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
            sigma  = LAM_i
          )[1]
          
          (max(rn_t, .SFA_CONSTANTS$MIN_POSITIVE) /
              max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
            exp(
              t(shift_vec) %*% ARR_i %*% e_i +
                0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
            )
        })
      },
      Ti       = t,
      ARR_i    = ARR,
      e_i      = e_i,
      LAM_i    = LAM,
      rd       = res_d,
      SIMPLIFY = FALSE
    )
    
    U <- unlist(U_list, use.names = FALSE)
    U <- pmin(U, 1)
    
    list(
      U = U,
      H = H,
      ridge_report = ridge_report
    )
  }
  
  
  ## ------------------------------------------------------------
  ## Finalize GTRE result object
  ##
  ## Handles:
  ##   - standard errors
  ##   - t-values
  ##   - result labeling
  ##   - sfareg return object
  ## ------------------------------------------------------------
  .gtre_finalize <- function(opt_obj,
                             prep,
                             data,
                             formula,
                             model_name,
                             call,
                             U,
                             H){
    
    # opt <- opt
    out <- out_template
    
    ## Standard errors
    if(isTRUE(any(opt$hessian == 0))){
      st_err <- rep(NA, length(opt$par))
    } else if(is.null(opt$hessian)) {
      st_err <- rep(NA, length(opt$par))
    } else {
      st_err <- suppressWarnings(
        sqrt(pmax(diag(solve(opt$hessian)), 0))
      )
    }
    
    t_val   <- opt$par / st_err
    out[1,] <- opt$par
    out[2,] <- st_err
    out[3,] <- t_val
    
    ## Rename intercept labels if needed
    if(length(colnames(out)[which(colnames(out) == "(Intercept)")]) > 2){
      colnames(out)[which(colnames(out) == "(Intercept)")] <-
        c("(Intercept x)", "(Intercept u)", "(Intercept h)")
    }
    
    total_time  <- End.Time
    
    results <- list(
      t(out),
      c(opt),
      data,
      total_time,
      start_v,
      model_name,
      formula,
      U,
      H,
      out["par",],
      out["st_err",],
      out["t-val",],
      call
    )
    
    class(results) <- "sfareg"
    names(results) <- c(
      "out",
      "opt",
      "data",
      "total_time",
      "start_v",
      "model_name",
      "formula",
      "U",
      "H",
      "coefficients",
      "std.errors",
      "t.values",
      "call"
    )
    
    results
  }
  
  
  ## ------------------------------------------------------------
  ## Main internal GTRE wrapper
  ##
  ## This is the single entry point for the GTRE branch.
  ##
  ## Recommended usage:
  ##   results <- .gtre(...)
  ##
  ## Then the higher-level package code can dispatch by model_name.
  ## ------------------------------------------------------------
  .gtre <- function(data,
                    individual,
                    y_var,
                    x_vars_vec,
                    z_vars,
                    zp_vars,
                    n_x_vars,
                    beta_hat,
                    beta_0,
                    beta_0_st,
                    sigma_v,
                    sigma_r,
                    start_val,
                    halton_num,
                    rand.gtre,
                    N,
                    inefdec_n,
                    maxit.bobyqa,
                    rand.psoptim,
                    maxit.psoptim,
                    PSopt,
                    maxit.optim,
                    optHessian,
                    Method,
                    formula,
                    call,
                    verbose = FALSE){
    
    # ## 2. Optimization
    # opt_obj <- .gtre_optimize(
    #   prep          = prep,
    #   inefdec_n     = inefdec_n,
    #   maxit.bobyqa  = maxit.bobyqa,
    #   rand.psoptim  = rand.psoptim,
    #   maxit.psoptim = maxit.psoptim,
    #   PSopt         = PSopt,
    #   maxit.optim   = maxit.optim,
    #   optHessian    = optHessian,
    #   Method        = Method,
    #   verbose       = verbose
    # )
    
    ## 3. Post-estimation technical efficiency
    te_obj <- .gtre_te(
      opt       = opt,
      prep      = prep,
      inefdec_n = inefdec_n
    )
    
    ## 4. Finalize result object
    .gtre_finalize(
      opt_obj    = opt_obj,
      prep       = prep,
      data       = data,
      formula    = formula,
      model_name = "GTRE_Z",
      call       = call,
      U          = te_obj$U,
      H          = te_obj$H
    )
  }
  
  results <- .gtre(
    data          = data,
    individual    = individual,
    y_var         = y_var,
    x_vars_vec    = x_vars_vec,
    z_vars        = z_vars,
    zp_vars       = zp_vars,
    n_x_vars      = n_x_vars,
    beta_hat      = beta_hat,
    beta_0        = beta_0,
    beta_0_st     = beta_0_st,
    sigma_v       = sigma_v,
    sigma_r       = sigma_r,
    start_val     = start_val,
    halton_num    = halton_num,
    rand.gtre     = rand.gtre,
    N             = N,
    inefdec_n     = inefdec_n,
    maxit.bobyqa  = maxit.bobyqa,
    rand.psoptim  = rand.psoptim,
    maxit.psoptim = maxit.psoptim,
    PSopt         = PSopt,
    maxit.optim   = maxit.optim,
    optHessian    = optHessian,
    Method        = Method,
    formula       = formula,
    call          = call,
    verbose       = verbose
  )
  return(results)
}
if(model_name ==     "TRE_Z")      {

delta          <- rep(0.1,length(z_vars))

if(isTRUE(is.numeric(start_val))){start_v <- start_val}
if(isFALSE(is.numeric(start_val))){start_v <- if(is.na(beta_0_st)){unname(c(sigma_v,sigma_r,beta_hat,delta))} else{unname(c(sigma_v,sigma_r,beta_0,beta_hat,delta))}}

out            <- matrix(0,nrow = 3,ncol = length(start_v))
rownames(out)  <- c("par","st_err","t-val") 
colnames(out)  <- c("sigv","sigr",colnames(data_x),z_vars)
      
if (isTRUE(is.numeric(start_val))){ start_v <- start_val}
if (isTRUE(is.numeric(halton_num))){R       <- halton_num}else{R <- ceiling(sqrt(nrow(data)))+100 }  ## Integral reps  

R_H     <- randtoolbox::halton(R+.SFA_CONSTANTS$HALTON_DISCARD,2,start = 1,normal = FALSE)[-c(1:.SFA_CONSTANTS$HALTON_DISCARD),c(1:2)]   
R_H     <- cbind( qnorm(R_H[,1]) , sqrt(2)* pracma::erfinv(R_H[,2]) )                ## using inverse error function for R_H2

if(!is.null(rand.gtre)){set.seed(rand.gtre)}   

mat <- matrix(0,nrow=R, ncol=9999)
for(v in 1:9999){mat[,v] <-  sample(R_H[,1])}

cor <- matrix(0,9999,1)
for(v in 1:9999){cor[v] <-  abs(cor(mat[,v],R_H[,2]))}

R_H <- cbind(mat[,which(cor==min(cor))] , R_H[,2])
rm(cor,v,mat)

# if(verbose){print(paste( "Primes 2 and 3 are in use, with 1,000 discards.  Correlation between R and H draws is:", round(cor(R_H)[1,2],10), sep = "" ),quote = FALSE) }
      
      indiv           <- noquote(as.vector(unique(data[,c(individual)])))
      t               <- rep(0, N)
      data_i <- Y <- eps <- data_i_vars <- data_z_vars <- R_h1 <- R_h2 <- as.list(rep(0,N))
      
      for (ii in 1:N) {
        data_i[[ii]]        <-  data[which(data[,c(individual)]==indiv[ii]),]
        t[ii]               <-  nrow(data_i[[ii]])
        R_h1[[ii]]          <-  t(matrix(rep(R_H[,1],t[[ii]]),R,t[[ii]]))
        R_h2[[ii]]          <-  abs(t(matrix(rep(R_H[,2],t[[ii]]),R,t[[ii]])))
        Y[[ii]]             <-  matrix(rep(data_i[[ii]][,y_var],R),t[[ii]],R)
        data_i_vars[[ii]]   <-  data.frame(data_i[[ii]][,c(x_vars_vec)] )
        data_z_vars[[ii]]   <-  data.frame(data_i[[ii]][,c(z_vars)])
        } 
      
fn <-  function(x){
x_x_vec    <- x[3:as.numeric(n_x_vars + 2)]

for (qq in 1:n_z_vars){
v              <- qq  + 2 + n_x_vars
z_z_vec[qq]    <- x[v]} 
  
fn1 = function(ii){ 
  eps[[ii]]     <- Y[[ii]]  - x[2]*R_h1[[ii]] 
  
  for (qq in 1:n_x_vars) {
  eps[[ii]] <- eps[[ii]] - x_x_vec[qq]*matrix(rep(data_i_vars[[ii]][,qq],R),t[[ii]],R)  }
  eps[[ii]] <- inefdec_n*eps[[ii]]
  
  sigma_u_fun  <- sqrt(exp(as.matrix(data_z_vars[[ii]])%*%z_z_vec))
  sigma_v_fun  <- x[1]
  sigma_fun    <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
  lamb_fun     <- sigma_u_fun/sigma_v_fun     
  
  prod_vec_n  <- log(mean(colProds((2/matrix(rep(sigma_fun,R),t[[ii]],R))* 
                 dnorm(eps[[ii]]/matrix(rep(sigma_fun,R),t[[ii]],R))*   
       pmax(pnorm(-eps[[ii]]*matrix(rep(lamb_fun,R),t[[ii]],R)/matrix(rep(sigma_fun,R),t[[ii]],R)),eps[[ii]]*0+.SFA_CONSTANTS$MIN_POSITIVE))))
  
  return(-prod_vec_n)
}

fn1_apply  <- unlist(lapply(1:N, fn1))

fn1_apply[which(fn1_apply==Inf)]   <-  (.SFA_CONSTANTS$MAX_VALUE)^.1
fn1_apply[which(fn1_apply==-Inf)]  <- -(.SFA_CONSTANTS$MAX_VALUE)^.1

return( sum( fn1_apply[is.finite(fn1_apply)] ) )} 

Start.Time <- start.time()

Opt.Bobyqa  <- opt.bobyqa(fn=fn, start_v=start_v, lower.bobyqa=-Inf, maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose)      
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 

differ<- 10
lower1    <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2) , start_v[-c(1:2)] -differ)

Opt.Psoptim <- opt.psoptim(fn=fn, start_v, lower.psoptim=lower1, upper.psoptim=c(start_v+differ), 
                           rand.psoptim=rand.psoptim,
                           maxit.psoptim=maxit.psoptim, psopt.TF=PSopt,verbose = verbose)
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

differ  <- 1
lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2) , start_v[-c(1:2)] -differ)

Opt.Optim <- opt.optim(fn = fn, start_v = start_v, lower.optim =lower1 ,upper.optim=c(start_v+differ), 
          maxit.optim=maxit.optim, opt.TF = optHessian, method=Method, optHessian= optHessian,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)         

if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

st_err     <- if (isTRUE(any(opt$hessian==0) | optHessian ==FALSE ) ){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)),0)))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
if(length(colnames(out)[which(colnames(out) == "(Intercept)")]) == 2){
colnames(out)[which(colnames(out) == "(Intercept)")] <- c("(Intercept x)" , "(Intercept u)")  }
      
## TE Measurements       
        NX    <- n_x_vars + 2
        NZ1   <- n_x_vars + 3
        NZ2   <- n_x_vars + n_z_vars + 2
        
        beta   <- opt$par[c(3:NX)]
        delta  <- opt$par[c(NZ1:NZ2)]
        
        sig_v  <- opt$par[1]
        sig_u  <- sqrt(exp((as.matrix(data.frame(data[,c(z_vars)] )))%*%delta))
        lamb   <- sig_u/sig_v
        sig    <- sqrt(sig_u^2 + sig_v^2)
        
        eps_hat    <- pmin( inefdec_n*(data[,y_var] - rowSums(t(t(data.frame(data[,c(x_vars_vec)] ))*beta))), 5)
        
        sig_star   <- (sig_u*sig_v)  /sig
        inner      <- (lamb*eps_hat) /sig
        U          <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/ pmax(1-pnorm(inner),.SFA_CONSTANTS$MIN_POSITIVE) )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
        U          <- pmax(U, 0)
        U          <- pmin(U, 1)
      
## Results 
results <- list(t(out),c(opt),data,End.Time,start_v,model_name,formula, U, out["par",], out["st_err",], out["t-val",], call)
class(results)  <- "sfareg"
names(results)  <- c("out","opt","data","total_time","start_v","model_name","formula","U", "coefficients", "std.errors", "t.values","call")
return(results)}    
if(model_name ==     "TFE")        {
Start.Tfe <- start.tfe(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars) 

data_i       <- Start.Tfe$data_i
data_i_vars  <- Start.Tfe$data_i_vars
eps          <- Start.Tfe$eps
I_t          <- Start.Tfe$I_t
I_t1         <- Start.Tfe$I_t1
indiv        <- Start.Tfe$indiv
one_t        <- Start.Tfe$one_t
one_t1       <- Start.Tfe$one_t1
out          <- Start.Tfe$out
start_v      <- Start.Tfe$start_v
t            <- Start.Tfe$t
upper        <- Start.Tfe$upper
Y            <- Start.Tfe$Y
  
like.tfe  <-  function(x){ 
x_x_vec    <- x[3:as.numeric(n_x_vars + 2)]  
    
fn1 = function(i){ 
      
eps_t    <- Y[[i]]  
      
for (qq in 1:n_x_vars) {
eps_t  <- eps_t - x_x_vec[qq]*demean(as.numeric(data_i_vars[[i]][,qq])) }
eps_t  <- eps_t*inefdec_n
eps_t1 <- eps_t[1:t[[i]]-1]
      
if(gamma==FALSE)   {E_t      <- ((x[1]*x[1])/t[[i]])*one_t[[i]]%*%t(one_t[[i]])}
else               {E_t      <- ( (x[1]/(1-x[1])) /t[[i]] )*one_t[[i]]%*%t(one_t[[i]])}
E_t1     <- (1/t[[i]])*one_t1[[i]]%*%t(one_t1[[i]])
l1       <- dmnorm(x= eps_t1,                  mean = rep(0,t[[i]]-1), varcov = x[2]*x[2]*(I_t1[[i]] - E_t1 ))
if(gamma==FALSE)   {l2       <- pmnorm(x= -(x[1]/x[2])*eps_t,  mean = rep(0,t[[i]]),       varcov =  I_t[[i]] + E_t )}
else               {l2       <- pmnorm(x= -(sqrt(x[1]/(1-x[1]))/x[2])*eps_t,  mean = rep(0,t[[i]]),  varcov =  I_t[[i]] + E_t )}
 
prod_vec_n  <- log(l1) + log(l2[1])  ## Log likelihood   
return(-prod_vec_n)
}
fn1_apply  <- unlist(lapply(1:N, fn1))
return(sum( fn1_apply[is.finite(fn1_apply)] ))}  

Start.Time <- start.time()

differ <- 2
Opt.Bobyqa <-  opt.bobyqa(fn=like.tfe, start_v=start_v, lower.bobyqa=c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2),rep(-Inf,n_x_vars)), 
               maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose)
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 

Opt.Psoptim <- opt.psoptim(fn=like.tfe, start_v, lower.psoptim=c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2), start_v[-c(1:2)]-differ), 
                           rand.psoptim = rand.psoptim, 
            upper.psoptim=c(start_v+differ), maxit.psoptim=maxit.psoptim, psopt.TF=PSopt,verbose = verbose)
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

Opt.Optim <- opt.optim(fn = like.tfe, start_v = start_v, lower.optim =c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2),rep(-Inf,n_x_vars)),
          upper.optim=c(start_v+differ), maxit.optim=maxit.optim, opt.TF = optHessian, method=Method, optHessian= optHessian,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)        

if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

beta             <- opt$par[-c(1:2)]  ## estimate the r_i's and then e(exp(u)|eps)'s 
lamb             <- opt$par[1]
sig              <- opt$par[2]
sig_u            <- (lamb*sig) / sqrt(1+lamb^2)
sig_v            <- sig_u/lamb
Y_mean           <- rep(0,N)
X_mean           <- matrix(0,N,ncol = length(x_vars_vec))
colnames(X_mean) <- x_vars_vec
      
for (ii in 1:N) {
        data_i[[ii]]  <- data[which(data[,c(individual)]==indiv[ii]),]
        Y_mean[ii]    <- mean(as.numeric(data_i[[ii]][,y_var]))
        X_mean[ii,]   <- colMeans(data.frame(data_i[[ii]][,c(x_vars_vec)] ))
      }
      
      r_hat_m      <- Y_mean - rowSums(t(beta*t(X_mean))) + inefdec_n*sqrt(2/pi)*sig_u
      r_hat_m_exp  <- rep(0,sum(t))
      t_cum        <- c(cumsum(t))
      
      for(m in 1: length(t)){
        B  <- t_cum[m]
        A  <- B +1 - t[m] 
        r_hat_m_exp[A:B] <- rep(r_hat_m[m],t[m])}
      
eps_hat    <- pmin(inefdec_n*(data[,y_var] - rowSums(t(t(data[,c(x_vars_vec)])*beta))- r_hat_m_exp),data[,y_var]*0)
sig_star   <- sig_u*sig_v/sig
inner      <- (lamb*eps_hat)/sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/(1-pnorm(inner)))*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
      
st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE ){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)),0)))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
results <- list(   t(out),c(opt),End.Time,   start_v,   r_hat_m,  exp_u_hat , model_name,  formula,  data, out["par",], out["st_err",], out["t-val",],call)
class(results) <- "sfareg"
names(results) <- c("out","opt","total_time","start_v","r_hat_m","exp_u_hat","model_name","formula","data","coefficients", "std.errors", "t.values","call")
return(results)
}
if(model_name ==     "FD")         {
Start.Time <- start.time()
      
if (isTRUE(is.numeric(start_val))) {start_v <- start_val} else{
        
        beta_hat         <- plm_fd$coefficients[c(x_vars_vec)]
        epsilon_hat      <- plm_fd$residuals
        beta_se          <- as.data.frame(summary(plm_fd)[1])$coefficients.Std..Error
        sfa_eps        <- pcs_c(Y  = as.numeric(epsilon_hat))[[1]]$par
        
        exp_u          <- sfa_eps[3]
        sigma_v        <- sqrt(sfa_eps[2]^2/ (1 + sfa_eps[1]^2))    #coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
        sigma_u        <- sigma_v*sfa_eps[1]                        #coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]
        
        sigmaSq_u        <- sigma_u^2
        sigmaSq_v        <- sigma_v^2
        mu               <- 0.1
        delta            <- rep(0.1,length(z_vars))
        start_v          <- unname(c(sigmaSq_u,sigmaSq_v,mu,beta_hat,delta))}
      
      
      out              <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)    <- c("par","st_err","t-val") 
      colnames(out)    <- c("sig_u2","sig_v2","mu",x_vars_vec,z_vars) 
      
      
      indiv           <- noquote(as.vector(unique(data[,c(individual)])))
      t               <- rep(0, N)
      data_i          <- as.list(rep(0,N))
      Y               <- as.list(rep(0,N))
      eps             <- as.list(rep(0,N))
      data_i_vars     <- as.list(rep(0,N))
      SIGMA           <- as.list(rep(0,N))
      data_z_vars     <- as.list(rep(0,N))
      
      for (ii in 1:N) {
        data_i[[ii]]              <- data[which(data[,c(individual)]==indiv[ii]),]
        t[ii]                     <- nrow(data_i[[ii]])
        Y[[ii]]                   <- as.numeric(data_i[[ii]][,y_var])
        data_i_vars[[ii]]         <- data.frame(data_i[[ii]][,c(x_vars_vec)] )
        data_z_vars[[ii]]         <- data.frame(data_i[[ii]][,c(z_vars)] )
        
        if(t[ii] == 2) {SIGMA[[ii]] = 2}
        if(t[ii] == 3) {SIGMA[[ii]] = matrix(c(2,-1,-1,2),2,2)} 
        if(t[ii] > 3) {
          SIGMA[[ii]]               <- matrix(0, nrow = t[ii]-1, ncol = t[ii]-1)
          diag(SIGMA[[ii]])         <-  2
          diag(SIGMA[[ii]][-1,  ])  <- -1
          diag(SIGMA[[ii]][  ,-1])  <- -1}}
      
like.fd <- function(x){ 
        
        for (q in 1:n_x_vars){
          v             <- q + 3
          x_x_vec[q]    <- x[v]}
        
        for (qq in 1:n_z_vars){
          m <- v + qq
          z_z_vec[qq]    <- x[m] }
        
        fn1 = function(i){
          
          eps_t    <- diff(Y[[i]], lag = 1) 
          eps_h    <- rep(0, t[i])
          
          for (qq in 1:n_x_vars) {eps_t <- eps_t - x_x_vec[qq] * diff(as.numeric(data_i_vars[[i]][,qq]), lag = 1)}
          for (qq in 1:n_z_vars) {eps_h <- eps_h + z_z_vec[qq] * as.numeric(data_z_vars[[i]][,qq])  }
          
          eps_h     <- diff(exp(eps_h), lag = 1)    
          eps_t     <- eps_t*inefdec_n            ## not exactly sure if this is sufficient for cost
          SIG       <- x[2]*SIGMA[[i]]
          sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h  + (1/x[1]))^-1
          mu_star   <- (  (x[3]/x[1]) - t(eps_t) %*% qr.solve(SIG) %*% eps_h  )*sig_star2
          
          l1     <-  - 0.5 * (t[i]-1) * log(2*pi) 
          l2     <-  - 0.5 * log(t[i])
          l3     <-  - 0.5 * (t[i]-1) * log(x[2])
          l4     <-  - 0.5 * t(eps_t) %*% qr.solve(SIG) %*% eps_t
          l5     <-    0.5 * ( (mu_star^2/sig_star2) -  ((x[3]^2) / x[1]) )
          l6     <-    log( sqrt(sig_star2)* max(pnorm(mu_star/ sqrt(sig_star2) ),.SFA_CONSTANTS$MIN_POSITIVE) )
          l7     <-  - log( sqrt(x[1]) *     max(pnorm( x[3] /  sqrt(x[1])      ),.SFA_CONSTANTS$MIN_POSITIVE) )
          
          
          prod_vec_n  <- sum(l1, l2, l3, l4, l5, l6, l7)      
          return(-prod_vec_n)}
        
        fn1_apply  <- unlist(lapply(1:N, fn1))
        
        return( sum( fn1_apply[is.finite(fn1_apply)] ) )}

differ  <- 2            

Opt.Bobyqa <- opt.bobyqa(fn=like.fd,  start_v=start_v, lower.bobyqa=c(rep(0.000001,2) , rep(-Inf,n_x_vars+1), rep(-Inf,length(z_vars))), 
              maxit.bobyqa=maxit.bobyqa, bob.TF=TRUE,verbose = verbose)     
start_v     <- Opt.Bobyqa$start_v
start_feval <- Opt.Bobyqa$start_feval
bob1        <- Opt.Bobyqa$bob1 

Opt.Psoptim <- opt.psoptim(fn=like.fd, start_v=start_v, lower.psoptim=c(rep(.SFA_CONSTANTS$MIN_POSITIVE,2), start_v[-c(1:2)]- differ), 
                           rand.psoptim = rand.psoptim, 
            upper.psoptim=c(start_v+differ), maxit.psoptim, psopt.TF=PSopt,verbose = verbose)
start_v     <- Opt.Psoptim$start_v
start_feval <- Opt.Psoptim$start_feval
opt00       <- Opt.Psoptim$opt00

Opt.Optim <- opt.optim(fn = like.fd, start_v = start_v, lower.optim =c(rep(0.000001,2) , rep(-Inf,n_x_vars+1), rep(-Inf,length(z_vars))),
          upper.optim=c(start_v+differ),maxit.optim=maxit.optim, opt.TF = optHessian, method=Method, optHessian= optHessian,verbose = verbose)
start_v     <- Opt.Optim$start_v
start_feval <- Opt.Optim$start_feval
opt         <- Opt.Optim$opt

End.Time <- end.time(Start.Time)   

if(optHessian==FALSE & PSopt == FALSE){opt <- bob1
st_err  <- rep(NA,length(opt$par))}

if(optHessian==FALSE & PSopt == TRUE){opt <- opt00
st_err  <- rep(NA,length(opt$par))}

## TE Measurements: ui's
NX    <- n_x_vars + 3
NZ1   <- n_x_vars + 4
NZ2   <- n_x_vars + n_z_vars + 3
beta  <- opt$par[c(4:NX)]
delta <- opt$par[c(NZ1:NZ2)]
sigu2 <- opt$par[1]
sigv2 <- opt$par[2]
mu    <- opt$par[3]
u_hat <- rep(0, sum(t))
h_hat <- exp(as.matrix(data.frame(subset(data,select = z_vars)))%*%delta)

for (i in 1:N) {
for (tt in 1:t[i]) {
eps_t     <- diff(Y[[i]], lag = 1) 
eps_h     <- rep(0, t[i]-1)

for (qq in 1:n_x_vars) {
m         <- 3 + qq
eps_t     <- eps_t - opt$par[m] * diff(as.numeric(data_i_vars[[i]][,qq]), lag = 1)}

for (qq in 1:n_z_vars) {
m         <- 3 + length(n_x_vars) + qq 
eps_h     <- eps_h + opt$par[m] * diff(as.numeric(data_z_vars[[i]][,qq]), lag = 1)}

eps_t     <- eps_t*inefdec_n
SIG       <- sigv2*SIGMA[[i]]
sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h  + (1/sigu2))^-1
mu_star   <- (  (mu/sigu2) - t(eps_t) %*% qr.solve(SIG) %*% eps_h  )*sig_star2 

num        <- if(i>1) { cumsum(t)[i-1]  + tt } else {tt}
u_hat[num]   <- h_hat[num] * (mu_star + (sqrt(sig_star2)*dnorm(mu_star/sqrt(sig_star2)) / max(pnorm(mu_star/sqrt(sig_star2)),.SFA_CONSTANTS$MIN_POSITIVE)   ) ) 
}}
      
exp_u_hat  <- exp(-u_hat)
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)
      
st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE){ rep(NA,length(opt$par)) }   else{suppressWarnings(sqrt(diag(qr.solve(opt$hessian))))} 
t_val      <- opt$par/st_err 
out[1,]        <- opt$par
out[2,]        <- st_err
out[3,]        <- t_val
results <- list(t(out),c(opt),End.Time,start_v,model_name,formula, u_hat, h_hat, exp_u_hat, data, out["par",], out["st_err",], out["t-val",],call)
class(results) <- "sfareg"
names(results) <- c("out","opt","total_time","start_v","model_name","formula", "u_hat", "h_hat", "exp_u_hat", "data", "coefficients", "std.errors", "t.values","call")
return(results)}
if(model_name ==     "GTRE_SEQ1")  {
Start.Time <- start.time()
      
## Sequential Method  
      sfa_eps       <- sfa(epsilon_hat   ~1, ineffDecrease = inefdec_TF)
      sfa_alp       <- sfa(alpha_hat     ~1, ineffDecrease = inefdec_TF)
      exp_eta       <- sfa_alp$mleParam[c(1)]
      exp_u         <- sfa_eps$mleParam[c(1)]
      sigma_u       <- coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]
      sigma_v       <- coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
      sigma_r       <- coef(sfa_alp,extraPar=TRUE)[c("sigmaU")]
      sigma_h       <- coef(sfa_alp,extraPar=TRUE)[c("sigmaV")]
      beta_0        <- beta_0_st + exp_u +exp_eta
      Lambda        <- sigma_u/sigma_v
      Sigma         <- sqrt(sigma_u^2 + sigma_v^2)
      sigma_se_r    <- NA
      sigma_r_se    <- NA
      sigma_h_se    <- NA
      lambda_se_r   <- NA
      beta_0_se     <- NA 
      gamma_uv      <- coef(sfa_eps,extraPar=TRUE)[c("gamma")]
      sigmaSq_uv    <- coef(sfa_eps,extraPar=TRUE)[c("sigmaSq")]
      gamma_hr      <- coef(sfa_alp,extraPar=TRUE)[c("gamma")]
      sigmaSq_hr    <- coef(sfa_alp,extraPar=TRUE)[c("sigmaSq")]
      gamma_uv_se   <- summary(sfa_eps)[[25]][c("gamma"),2]
      sigmaSq_uv_se <- summary(sfa_eps)[[25]][c("sigmaSq"),2]
      gamma_hr_se   <- summary(sfa_alp)[[25]][c("gamma"),2]
      sigmaSq_hr_se <- summary(sfa_alp)[[25]][c("sigmaSq"),2]
      
      other_parms            <- as.matrix(c(Lambda,Sigma,beta_0))
      rownames(other_parms)  <- c("lambda","sigma","beta_0")
      start_v                <- unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_hat))
      out                    <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)          <- c("par","st_err","t-val") 
      colnames(out) <- if(isTRUE(intercept==0)) {c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )} else{c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x)[-c(1)] )}
      
     End.Time <- end.time(Start.Time)   
      st_err     <- rep(NA,ncol(out)) 
      st_err     <- if(isTRUE(intercept==0)) {c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[,2])} else{c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[-c(1),2])}
      t_val      <- start_v/st_err
      out[1,]    <- start_v
      out[2,]    <- st_err
      out[3,]    <- t_val
results <- list(   t(out),End.Time,     other_parms,  model_name,  formula,   data,   out["par",],    out["st_err",], out["t-val",],call)
class(results) <- "sfareg"
names(results) <- c("out","total_time","other_parms","model_name","formula", "data", "coefficients", "std.errors",   "t.values","call")
return(results)}
if(model_name ==     "GTRE_SEQ2")  {
Start.Time <- start.time()
## Sequential Method following 1995 paper  
## take second and third moments of alpha_hat and epsilon_hat
      alp_2m <- mean(alpha_hat^2)
      alp_3m <- min(0,mean(alpha_hat^3))
      eps_2m <- mean(epsilon_hat^2)
      eps_3m <- min(0,mean(epsilon_hat^3))
      
      gamma_uv   <- min(1, 1/(eps_2m*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3) + (2/pi)  ))
      gamma_hr   <- min(1, 1/(alp_2m*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3) + (2/pi)  ))
      sigmaSq_uv <- eps_2m + (2/pi)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(2/3)
      sigmaSq_hr <- alp_2m + (2/pi)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(2/3)
      beta_0     <- beta_0_st + sqrt(2/pi)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(1/3) + sqrt(2/pi)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(1/3)
      
      ## calculate the ten needed central moments 
      mu_2_eps   <- sigmaSq_uv*((1-gamma_uv) + gamma_uv*((pi -2)/pi) )
      mu_3_eps   <- sigmaSq_uv^(3/2)*(sqrt(2/pi)*(1-(4/pi))*gamma_uv^(3/2))
      mu_4_eps   <- sigmaSq_uv^2*(3*(1-gamma_uv)^2  + ((6*(pi-2)*gamma_uv*(1-gamma_uv))/pi) + gamma_uv^2*(3- (4/pi) - (12/pi^2)) )
      mu_5_eps   <- sigmaSq_uv^(5/2)*gamma_uv^(3/2)*sqrt(2/pi)*(10*(1-(4/pi))*(1-gamma_uv) + (7-(20/pi)-(16/pi^2))*gamma_uv )
      mu_6_eps   <- sigmaSq_uv^3*( 15*(1-gamma_uv)^3   +  
                                     (45*(pi-2)*(1-gamma_uv)^2*gamma_uv/pi)  +
                                     15*(3 - (4/pi) - (12/pi^2))*(1-gamma_uv)*gamma_uv^2  +
                                     (15  -  (6/pi) - (100/pi^2)  - (40/pi^3)  )*gamma_uv^3         )
      
      mu_2_alp   <- sigmaSq_hr*((1-gamma_hr) + gamma_hr*((pi -2)/pi) )
      mu_3_alp   <- sigmaSq_hr^(3/2)*(sqrt(2/pi)*(1-(4/pi))*gamma_hr^(3/2))
      mu_4_alp   <- sigmaSq_hr^2*(3*(1-gamma_hr)^2  + ((6*(pi-2)*gamma_hr*(1-gamma_hr))/pi) + gamma_hr^2*(3- (4/pi) - (12/pi^2)) )
      mu_5_alp   <- sigmaSq_hr^(5/2)*gamma_hr^(3/2)*sqrt(2/pi)*(10*(1-(4/pi))*(1-gamma_hr) + (7-(20/pi)-(16/pi^2))*gamma_hr )
      mu_6_alp   <- sigmaSq_hr^3*( 15*(1-gamma_hr)^3   +  
                                     (45*(pi-2)*(1-gamma_hr)^2*gamma_hr/pi)  +
                                     15*(3 - (4/pi) - (12/pi^2))*(1-gamma_hr)*gamma_hr^2  +
                                     (15  -  (6/pi) - (100/pi^2)  - (40/pi^3)  )*gamma_hr^3         )
      
      var_2m_eps  <-  (1/nrow(data))*(mu_4_eps - mu_2_eps^2)
      var_3m_eps  <-  (1/nrow(data))*(mu_6_eps - mu_3_eps^3 -6*mu_2_eps*mu_4_eps + 9*mu_2_eps) 
      cov_23m_eps <-  (1/nrow(data))*(mu_5_eps - 4*mu_2_eps*mu_3_eps)
      
      var_2m_alp  <-  (1/N)*(mu_4_alp - mu_2_alp^2)
      var_3m_alp  <-  (1/N)*(mu_6_alp - mu_3_alp^3 -6*mu_2_alp*mu_4_alp + 9*mu_2_alp) 
      cov_23m_alp <-  (1/N)*(mu_5_alp - 4*mu_2_alp*mu_3_alp)
      
      ## define 8 needed derivatives 
      d_beta_0_d_m3_eps  <- (pi/(pi-4))*(1/3)    *       (sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3)
      d_sigma2_d_m3_eps  <- sqrt(2/pi)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-1/3)
      d_gamma_d_m2_eps   <- -(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3)*(eps_3m*(  sqrt(pi/2)*(pi/(pi-4))*eps_3m    )  + (2/pi))
      d_gamma_d_m3_eps   <- eps_2m*sqrt(pi/2)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-5/3)*(eps_3m*(sqrt(pi/2)*(pi/(pi-4))*eps_3m )^(-2/3) + (2/pi) )^(-2)
      
      d_beta_0_d_m3_alp  <- (pi/(pi-4))*(1/3)    *       (sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3)
      d_sigma2_d_m3_alp  <- sqrt(2/pi)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-1/3)
      d_gamma_d_m2_alp   <- -(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3)*(alp_3m*(  sqrt(pi/2)*(pi/(pi-4))*alp_3m    )  + (2/pi))
      d_gamma_d_m3_alp   <- alp_2m*sqrt(pi/2)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-5/3)*(alp_3m*(sqrt(pi/2)*(pi/(pi-4))*alp_3m )^(-2/3) + (2/pi) )^(-2)
      
      beta_0_se     <- sqrt(beta_se[1]  + d_beta_0_d_m3_eps^2*var_3m_eps + d_beta_0_d_m3_alp^2*var_3m_alp )
      
      sigmaSq_uv_se <- sqrt(var_2m_eps + d_sigma2_d_m3_eps^2*var_3m_eps  + d_sigma2_d_m3_eps * cov_23m_eps)
      gamma_uv_se   <- sqrt(d_gamma_d_m2_eps^2*var_2m_eps  + d_gamma_d_m3_eps^2* var_3m_eps + d_gamma_d_m2_eps*d_gamma_d_m3_eps*cov_23m_eps     )
      
      sigmaSq_hr_se <- sqrt(var_2m_alp + d_sigma2_d_m3_alp^2*var_3m_alp  + d_sigma2_d_m3_alp * cov_23m_alp)
      gamma_hr_se   <- sqrt(d_gamma_d_m2_alp^2*var_2m_alp  + d_gamma_d_m3_alp^2* var_3m_alp + d_gamma_d_m2_alp*d_gamma_d_m3_alp*cov_23m_alp     )
      
      start_v        <- if(is.na(beta_0_st)) {unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_hat))} else {unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_0,beta_hat))}
     End.Time <- end.time(Start.Time)   
      out            <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)  <- c("par","st_err","t-val") 
      colnames(out)  <- if(isTRUE(intercept==0)) {c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )} else{c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )}
      
      st_err     <- rep(NA,ncol(out)) 
      st_err     <- if(is.na(beta_0_st)){c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[c(x_vars_vec),2])} else {c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,beta_0_se,summary(plm_gtre)$coefficients[-c(1),2])}
      t_val      <- start_v/st_err
      out[1,]    <- start_v
      out[2,]    <- st_err
      out[3,]    <- t_val

      ## look at individual sigmas 
      Sigma_u <- sqrt(gamma_uv*sigmaSq_uv)    
      Sigma_v <- sqrt((1-gamma_uv)*sigmaSq_uv) 
      Sigma_h <- sqrt(gamma_hr*sigmaSq_hr)   
      Sigma_r <- sqrt((1-gamma_hr)*sigmaSq_hr) 
      Lambda  <- sigma_u/sigma_v
      Sigma   <- sqrt(sigmaSq_uv)
      
      other_parms    <- as.matrix(c(Sigma_u,Sigma_v,Sigma_h, Sigma_r, Lambda, Sigma))
      rownames(other_parms) <- c("sigma_u","sigma_v","sigma_h", "sigma_r","lambda","sigma")
      
      results <- list(   t(out),End.Time,     other_parms,  model_name,  formula,   data,   out["par",],    out["st_err",], out["t-val",],call)
      class(results) <- "sfareg"
      names(results) <- c("out","total_time","other_parms","model_name","formula", "data", "coefficients", "std.errors",    "t.values","call")
      return(results)}

else {return(c("This is not a valid command"))}}
