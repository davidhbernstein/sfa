## ------------------------------------------------------------
## S3 methods for class "sfareg"
##
## psfm(), sfm(), zsfm(), and ttsfm() all return objects of class
## "sfareg", but previously only print.sfareg() and summary.sfareg()
## existed (see data.processing.R). Users had to reach into
## object$coefficients / object$opt by hand for anything else. These
## add the standard R modeling generics so sfareg objects behave like
## other model objects (coef(fit), logLik(fit), AIC(fit), BIC(fit),
## vcov(fit)).
##
## Sign convention: every estimator in this package minimizes its
## objective function (bobyqa/psoptim/optim all minimize by default;
## see opts.R, which does not flip the sign of the user-supplied fn),
## and each fn is written to return the NEGATIVE summed
## log-likelihood -- so logLik = -object$opt$value. This matches
## print.sfareg()/summary.sfareg(), which already display
## -object$opt$value as "log likelihood".
##
## psfm()'s GTRE_SEQ1/GTRE_SEQ2 branches are moment-based (not MLE)
## and carry no $opt component at all; logLik/AIC/BIC are undefined
## for those fits and will return NA with a warning rather than
## erroring.
## ------------------------------------------------------------

coef.sfareg <- function(object, ...){
  object$coefficients
}

vcov.sfareg <- function(object, ...){
  p  <- length(object$coefficients)
  nm <- names(object$coefficients)

  if(!is.null(object$opt) && !is.null(object$opt$hessian)){
    V <- tryCatch(solve(object$opt$hessian), error = function(e) NULL)
    if(!is.null(V)){
      dimnames(V) <- list(nm, nm)
      return(V)
    }
    warning("Hessian stored on this fit could not be inverted; falling back to a diagonal approximation built from the reported standard errors.", call. = FALSE)
  }

  if(!is.null(object$std.errors) && !all(is.na(object$std.errors))){
    V <- diag(object$std.errors^2, nrow = p)
    dimnames(V) <- list(nm, nm)
    return(V)
  }

  warning("No Hessian or standard errors are available on this fit (was optHessian = FALSE?); returning a matrix of NAs.", call. = FALSE)
  matrix(NA_real_, p, p, dimnames = list(nm, nm))
}

logLik.sfareg <- function(object, ...){
  if(is.null(object$opt) || is.null(object$opt$value)){
    warning("This fit has no stored optimizer output (e.g. psfm()'s GTRE_SEQ1/GTRE_SEQ2 are moment-based, not maximum likelihood), so logLik() is not defined for it.", call. = FALSE)
    return(NA_real_)
  }
  val <- -object$opt$value
  attr(val, "df")   <- length(object$coefficients)
  attr(val, "nobs") <- nobs.sfareg(object)
  class(val) <- "logLik"
  val
}

nobs.sfareg <- function(object, ...){
  if(!is.null(object$data)){
    return(nrow(object$data))
  }
  tryCatch({
    dcall <- object$call$data
    if(is.null(dcall)) return(NA_integer_)
    nrow(eval(dcall, envir = parent.frame()))
  }, error = function(e) NA_integer_)
}
