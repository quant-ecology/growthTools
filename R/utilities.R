

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

#' Extract number of observations from mle2 object
#' 
#' By default, mle2 objects from bbmle produced by fitting a model to data without
#' an explicit formula (e.g., using a black-boxed NLL calculator function instead)
#' do not provide the number of observations underlying the fit. This extends the 
#' nobs() generic with a specific method for mle2 objects.
#' 
#' @export
nobs.mle2<-function(object){
  length(object@data[[1]])
}


#' Make sure Julia is available if calling a function that needs it
#' 
ensure_julia <- function() {
  
  # need JuliaCall package
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("JuliaCall not installed")
  }
  
  # detect Julia session
  julia_alive <- !inherits(
    try(JuliaCall::julia_eval("1"), silent = TRUE),
    "try-error"
  )
  if (!julia_alive) {
    JuliaCall::julia_setup(installJulia = FALSE)
  }
  
  # check if already fully initialized (R-side AND Julia-side)
  julia_ready <- tryCatch({
    JuliaCall::julia_eval("isdefined(Main, :solve)")
  }, error = function(e) FALSE)
  if (exists(".julia_initialized", envir = .GlobalEnv, inherits = FALSE) &&
      julia_ready) {
    return(invisible(NULL))
  }
  
  # otherwise, julia is set up but has not yet been initialized for ODE work, so:      
  # HARD LOAD SciML STACK (critical)
  JuliaCall::julia_eval("using Pkg; using DifferentialEquations;")
  
  # VERIFY it actually worked (important)
  JuliaCall::julia_eval("ODEProblem")
  
  # ODE in log-space
  JuliaCall::julia_eval("
    function f!(du,u,p,t)
      alpha, vmax, c, d = p

      r = u[1]
      n = u[2]   # log(N)

      du[1] = alpha * exp(n) * (c*d - r/(r+1))
      du[2] = vmax * (r/(r+1) - d)
    end
  ")
  
  # evaluate a dummy/template version of the problem:
  #JuliaCall::julia_eval("u0 = [1.0, 1.0]; p  = [1.0, 1.0, 0.5, 0.5]; tspan = (0.0, 10.0); prob_template = ODEProblem(f!, u0, tspan, p)")
  
  assign(".julia_initialized", TRUE, envir = .GlobalEnv)
}


#' Define global parameters
r0.global <<- 100

