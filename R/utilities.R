

#' Calculate R2
#' 
#' 1-sum((obs-pds)^2)/sum((obs-mean(obs))^2)
#' 
#' @param pds Predicted y values (could arise from predict(model) for an mle2 model)
#' @param obs Observed y values
#' 
#' @return R2 value
#' 
#' @export
get.R2<-function(pds,obs){
  #pds<-predict(model)
  #cor(pred, obs)^2
  rsqr<-1-sum((obs-pds)^2)/sum((obs-mean(obs))^2)
  rsqr
}


#' Delta method function (modified from emdbook)
#' 
#' See `?deltavar()` in package emdbook. Only change was to replace `expr <- as.expression(substitute(fun))` with `expr<-fun`. This was necessary b/c of problems that arose when
#' invoking `deltavar()` inside of another function. Specifically, I wanted to call 
#' `deltavar()` each time an outer function was executed, but apply deltavar to a unique
#' set of x-values upon each execution of the outer function. For reasons I don't 
#' entirely understand, this range of x-values defined within the outer function's 
#' environment was not being adequately passed to the internal environment of deltavar.
#' There's probably a classier/more resilient way of fixing this, but in the short 
#' term, I patched it over by creating `deltavar2()`, and using parse/a custom string in
#' `get.nbcurve.tpc()` when `deltavar2()` is called.
#' 
#' @param fun Function to calculate the variance of, given parameter estimates in meanvals
#' @param vars 'list of variable names: needed if params does not have names, or if some of the values specified in params should be treated as constant'
#' @param meanval 'possibly named vector of mean values of parameters'
#' @param Sigma 'numeric vector of variances or variance-covariance matrix'
#' @param verbose 'print details?'
#' 
#' @export
deltavar2<-function (fun, meanval = NULL, vars, Sigma, verbose = FALSE) 
{
  #expr <- as.expression(substitute(fun))    
  expr<-fun  # Changed by CTK
  nvals <- length(eval(expr, envir = as.list(meanval)))
  vecexp <- nvals > 1
  if (missing(vars)) {
    if (missing(meanval) || is.null(names(meanval))) 
      stop("must specify either variable names or named values for means")
    vars <- names(meanval)
  }
  derivs <- try(lapply(vars, stats::D, expr = expr), silent = TRUE)
  symbderivs <- TRUE
  if (inherits(derivs, "try-error")) {
    symbderivs <- FALSE
    warning("some symbols not in derivative table, using numeric derivatives")
    nderivs <- with(as.list(meanval), numericDeriv(expr[[1]], 
                                                   theta = vars))
    nderivs <- attr(nderivs, "gradient")
  }
  else {
    nderivs <- sapply(derivs, eval, envir = as.list(meanval))
  }
  if (verbose) {
    if (symbderivs) {
      cat("symbolic derivs:\n")
      print(derivs)
    }
    cat("value of derivs:\n")
    print(nderivs)
  }
  if (!is.matrix(Sigma) && length(Sigma) > 1) 
    Sigma <- diag(Sigma)
  if (vecexp && is.list(nderivs)) 
    nderivs <- do.call("cbind", nderivs)
  if (is.matrix(nderivs)) { 
    # turning off this transformation of nderivs, which was causing predict.monod to fail
    # caution: test this for predict.tpc, which motivated this line in the first place...
    #nderivs<-t(matrix(nderivs[,(names(meanval) %in% colnames(Sigma))])) # added by CTK
    r <- apply(nderivs, 1, function(z) c(z %*% Sigma %*% 
                                           matrix(z)))
  }
  else r <- c(nderivs %*% Sigma %*% matrix(nderivs))
  r
}

#' Approximate 2nd derivative (finite difference)
#' 
#' @param fx A vector of 3 function values, f(x-h), f(x), f(x+h)
#' @param h The step size used in the finite difference calculation
#' 
#' @export
fd2.central<-function(fx,h){
  (fx[3]-2*fx[2]+fx[1])/(h^2)
}


#' Start assuming Julia is not initialized
.julia_initialized <- FALSE

#' Make sure Julia is available if calling a function that needs it
#' 
ensure_julia <- function() {
  
  if (!.julia_initialized) {
    JuliaCall::julia_setup(installJulia = FALSE, quiet = TRUE)
    .julia_initialized <<- TRUE
  }
}

