#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Growth rate estimation routines:   ####

# Developed by CTK for NSF Dimensions project

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#' Helper function for smoothed lagged/saturating abundance equations
#' 
#' @param x time
#' @param b slope
#' @param s smoothing parameter
sqfunc<-function(x,b,s){
  (1/2)*sqrt(b*(4*s+b*x^2))
}


#' Equations for modeling abundance time series
#' 
#' Intended to allow the extraction of exponential growth rates from time series 
#' while accounting for the presence of initial lags in growth, saturating abundances,
#' or both in the same time series. These equations provide smoothed piecewise linear
#' functions, where lagged or saturated portions of the time series maintain constant 
#' abundance, and elsewhere log(abundance) increases linearly.
#' 
#' This approach for lag, saturation, and lag+saturation are based on:
#' https://stats.stackexchange.com/questions/149627/piecewise-regression-with-constraints
#' which invokes a smooth approximation to a piecewise linear function,
#' where parameter s determines the smoothness around break-points. Generally, as s->0 
#' this smooth model approximates more closely the piecewise linear one. The s term 
#' could be fit explicitly, but for now it is fixed at a small number (1E-10).
#' 
#' Note: Currently, only the linear model without lag or saturation can produce negative
#' growth rate estimates. The lagged/saturating models will be extended to allow this 
#' possibility in future versions of this package.
#' 
#' @param x Time variable
#' @param a Initial log(abundance) at time = 0
#' @param b slope of the increasing linear portion of the time series, must be >=0
#' @param b2 slope of the decreasing linear portion of the time series (satdecay), must be <=0
#' @param B1 Time point where abundance starts to increase (leaves lag phase)
#' @param B2 Time point where abundance stops increasing (saturates)
#' @param s Smoothing parameter; as this term -> 0, these continuous functions approach true piecewise equations
#' 
#' @return log(abundance) at time x as a function of model parameters
#' 
#' @examples 
#' 
#' curve(lag(x,5,1,4,s=1E-10),0,10,col='green',ylim=c(0,11),ylab='Abundance')
#' curve(sat(x,0.9,1,8,s=1E-10),0,10,col='red',add=TRUE)
#' curve(lagsat(x,5.1,1,4,8,s=1E-10),0,10,col='blue',add=TRUE)
#' curve(flr(x,10,-1,8,s=1E-10),0,10,col='purple',add=TRUE)
#' 
#' @export
lagsat<-function(x,a,b,B1,B2,s=1E-10){
  a + (1/2)*b*(B2-B1) + sqfunc(B1-x,b,s) - sqfunc(B2-x,b,s)
}

#' @describeIn lagsat Lagged increasing linear function
#' @export
lag<-function(x,a,b,B1,s=1E-10){
  sqfunc(B1-x,b,s)-(b/2)*(B1-x)+a
}

#' @describeIn lagsat Saturating linear function
#' @export
sat<-function(x,a,b,B2,s=1E-10){
  a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s)
}

#' @describeIn lagsat Saturating then decaying linear functions
#' @export
satdecay<-function(x,a,b,b2,B2,s=1E-10){
  a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s) - (-b2/2)*(x-B2) - sqfunc(B2-x,-b2,s)
}

#' @describeIn lagsat Floored decreasing linear function
#' @export
flr<-function(x,a,b,B2,s=1E-10){
  b <- -1*b
  a - (1/2)*(b)*(B2) - sqfunc(-x,b,s) + sqfunc(B2-x,b,s)
}

#' ODE Equations for modeling abundance time series
#' 
#' Intended to allow the extraction of exponential growth rates from time series 
#' while accounting for unobserved, resource-dependent growth, saturating abundance,
#' and potentially declines in abundance after saturation due to accumulation of
#' resource in a recalcitrant pool. The satdecay() function is a piecewise linear
#' approximation of these dynamics.
#' 
#' NOTE: currently, this function depends on computation of numerical ODE solution
#' in Julia (for computational efficiency), and will not work without a functional
#' connection to Julia. 
#' 
#' @param x Time variable
#' @param alpha Affinity; defined as vmax over k, must be >0
#' @param vmax Maximum uptake rate, must be >0
#' @param cpar Proportion of dead biomass returned to the labile resource pool, must be 0 <= cpar <= 1
#' @param dpar Mortality rate, must be >0
#' @param r0 Initial resource concentration at time t=0 (unitless, scaled by k, must be >0)
#' @param n0 Initial log(abundance) at time t=0
#' 
#' @return log(abundance) at time x as a function of model parameters
#' 
#' @export
satdecay.ode <- function(x, alpha, vmax, cpar, dpar, r0, n0) {
  
  # make sure Julia is accessible
  ensure_julia()
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("JuliaCall is required for fitting satdecay_ode model")
  }
  
  # define time range
  tmax <- max(max(x), 10)
  
  julia_assign("times", x)
  julia_assign("p_new", c(alpha, vmax, cpar, dpar))
  julia_assign("u0_new", c(r0, n0))
  julia_assign("tmax", tmax)
  
  # below only works if x is more than one value
  vals <- julia_eval("prob = remake(prob_template,u0=u0_new,p=p_new,tspan=(0.0, tmax)); sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6); [sol(t)[2] for t in times]")
  return(vals)
}

#' Estimate timing of peak abundance (satdecay ODE model)
#' 
#' Given a set of estimated coefficients for a satdecay ODE model, use numerical
#' methods in Julia to calculate the timing of peak ln(abundance).
#' 
#' @param cfs Coefficients of satdecay ODE model, e.g. from get.gr.satdecay.ode()
#' @param r0 Initial resource concentration (set arbitrarily to 10 throughout)
#' @param tmax Maximum of time domain
#' 
#' @export
satdecay.ode.peak.time <- function(cfs, r0=10, tmax=100){
  # make sure Julia is accessible
  ensure_julia()
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("JuliaCall is required for fitting satdecay_ode model")
  }
  
  julia_assign("p_new", c(cfs$alpha, cfs$vmax, cfs$c, cfs$d))
  julia_assign("u0_new", c(r0, cfs$n0))
  julia_assign("tmax_local", tmax)
  
  julia_eval("prob = remake(prob_template,u0=u0_new,p=p_new,tspan=(0.0,tmax_local)); sol = solve(prob,Tsit5(),saveat=0.01,reltol=1e-8,abstol=1e-8,save_everystep=false); rvals = [u[1] for u in sol.u]; thresh = p_new[4] / (1 - p_new[4]); idx = findfirst(x -> x <= thresh, rvals); idx === nothing ? NaN : sol.t[idx];")
}

#' Extract exponential growth rate assuming exponential growth
#' 
#' This function fits a linear model to ln(abundance) data.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a linear model regressing ln(abundance) on time
#' 
#' @export
get.gr.linear<-function(x,y,plotQ=F,fpath=NA,id=''){
  lm1<-stats::lm(y~x)
  
  if(plotQ){
    if(!is.na(fpath)){
      grDevices::pdf(fpath)
      graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      graphics::abline(lm1,col='red')
      grDevices::dev.off()
    }else{
      graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      graphics::abline(lm1,col='red')
    }
  }
  
  return(lm1)
}

#' Extract exponential growth rate assuming lagged exponential growth
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances are nearly constant for several time points, before 
#' exponential growth kicks in.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
get.gr.lag<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  #slopes <- zoo::rollapply(data, 3, localslope, by.column=F)
  
  fit.lag<-try(nlsLM(y ~ lag(x,a,b,B1,s=1E-10),
                 start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), b=1),data = data,
                 lower = c(B1=-Inf,a=-Inf,b=0.0001),
                 control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  if(class(fit.lag)=='try-error'){
    fit.lag<-try(nlsLM(y ~ lag(x,a,b,B1,s=1E-10),
                       start = c(B1=10, a=min(y), b=1),data = data,
                       lower = c(B1=-Inf,a=-Inf,b=0.0001),
                       control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  }
  if(class(fit.lag)=='try-error'){
    if(!grepl(attr(fit.lag,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.lag,"condition"))
    }
    #print('fit.lag failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.lag)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),min(x),max(x),n = 400,add=TRUE,col='blue')
        graphics::curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),cfs$B1,max(x),n = 400,add=TRUE,col='red')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),min(x),max(x),n = 400,add=TRUE,col='blue')
        graphics::curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),cfs$B1,max(x),n = 400,add=TRUE,col='red')
      }
    }
  }

  return(fit.lag)
}

#' Extract exponential growth rate assuming exponential growth that saturates
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances increase linearly at first, but then saturate and
#' remain constant.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
get.gr.sat<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  slopes <- zoo::rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  a.guess<-coef(stats::lm(y~x))[[1]]
  
  fit.sat<-try(nlsLM(y ~ sat(x,a,b,B2,s=1E-10),
                 start=c(B2=mean(x)+(max(x)-mean(x))/2,a=a.guess,b=round(max(c(slopes,0.0001)),5)),
                 data = data,
                 lower = c(B2=-Inf,a=-Inf,b=0.0001),
                 control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  if(class(fit.sat)=='try-error'){
    fit.sat<-try(nlsLM(y ~ sat(x,a,b,B2,s=1E-10),
                       start=c(B2=10,a=a.guess,b=round(max(c(slopes,0.0001)),5)),
                       data = data,
                       lower = c(B2=-Inf,a=-Inf,b=0.0001),
                       control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  }
  if(class(fit.sat)=='try-error'){
    if(!grepl(attr(fit.sat,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.sat,"condition"))
    }
    #print('fit.sat failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.sat)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
      }
    }
  }  
  
  return(fit.sat)
}


#' Extract exponential growth rate assuming exponential growth that saturates and 
#' then decays
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances increase linearly at first, but then saturate and
#' then abruptly decay linearly instead of remaining constant.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
get.gr.satdecay<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  slopes <- zoo::rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  
  #a.guess <- coef(stats::lm(y ~ x))[[1]]
  a.guess <- mean(y[1:2])
  
  #B2.guess<-mean(x) + (max(x) - mean(x))/2
  B2.guess<-x[which(y==max(y)[1])]
  
  fit.satdecay<-try(nlsLM(y ~ satdecay(x,a,b,b2,B2,s=1E-10),
                     start=c(B2=B2.guess,a=a.guess,
                             b=round(max(c(slopes,0.0001)),5),
                             b2=round(min(c(slopes,-0.0001)),5)),
                     data = data,
                     lower = c(B2=-Inf,a=-Inf,b=0.0001,b2=-Inf), # do we want to constrain b2?
                     control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  if(class(fit.satdecay)=='try-error'){
    fit.satdecay<-try(nlsLM(y ~ satdecay(x,a,b,b2,B2,s=1E-10),
                       start=c(B2=10,a=a.guess,b=round(max(c(slopes,0.0001)),5),
                               b2=round(min(c(slopes,-0.0001)),5)),
                       data = data,
                       lower = c(B2=-Inf,a=-Inf,b=0.0001,b2=-Inf),
                       control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  }
  if(class(fit.satdecay)=='try-error'){
    if(!grepl(attr(fit.satdecay,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.satdecay,"condition"))
    }
    #print('fit.satdecay failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.satdecay)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay(x,cfs$a,cfs$b,cfs$b2,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(satdecay(x,cfs$a,cfs$b,cfs$b2,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay(x,cfs$a,cfs$b,cfs$b2,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(satdecay(x,cfs$a,cfs$b,cfs$b2,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
      }
    }
  }  
  
  return(fit.satdecay)
}


#' Extract exponential growth rate assuming exponential growth that saturates and 
#' then decays using full ODE fitting
#' 
#' This function fits ODE model to ln(abundance) data, with 
#' the assumption that abundances increase linearly at first, as resource is 
#' depleted in the background. If resources are perfectly recycled into the 
#' resource pool, abundances will saturate. But, imperfect recycling leads to
#' subsequent decay in abundance as cells die, and recalcitrant nutrients 
#' accumulate. Under particular circumstances, this dynamic is well approximated
#' by a piecewise linear function; see get.gr.satdecay(). Currently, this 
#' functionality only works if R is connected to Julia, which provides fast
#' ODE solutions
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns an mle2 regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
get.gr.satdecay.ode<-function(x,y,plotQ=F,fpath=NA,id=''){
  data<-data.frame(x=x,y=y)
  
  if(length(unique(x))<4){
    stop("error: fewer than four distinct time steps provided to get.gr.satdecay.ode")
  }
  
  ## Formulate starting guesses:
  
  # initial slope should be ~= vmax*(1-d)
  fit_early <- lm(y[1:3] ~ x[1:3]) 
  slope0 <- coef(fit_early)[2]
  
  # final slope should be ~= -vmax*d
  fit_late <- lm(y[(length(x)-2):length(x)] ~ x[(length(x)-2):length(x)])
  slope_end <- coef(fit_late)[2]
  
  vmax.guess <- slope0 - slope_end
  d.guess    <- -slope_end / vmax.guess
  n0.guess <- y[which.min(x)]
  c.guess <- 0.2
  alpha.guess <- 0.1 * vmax.guess # careful with this one; linked to r0 assumption
  
  # set up for likelihood calculation:
  JuliaCall::julia_assign("times_obs", x)
  JuliaCall::julia_assign("tmax_global", max(max(x), 10))

  # precompile solver: (is this necessary/helpful?)
  JuliaCall::julia_eval("prob = remake(prob_template); sol = solve(prob, Tsit5(), saveat=times_obs); nothing")
  
  # local version of satdecay.ode(), to optimize run time. Uses fixed time vals
  satdecay.ode.local <- function(x, alpha, vmax, cpar, dpar, r0, n0) {
    # define time range
    tmax <- max(max(x), 10)
    
    JuliaCall::julia_assign("p_new", c(alpha, vmax, cpar, dpar))
    JuliaCall::julia_assign("u0_new", c(r0, n0))
    
    # below only works if x is more than one value
    vals <- JuliaCall::julia_eval("prob = remake(prob_template,u0=u0_new,p=p_new,tspan=(0.0, tmax_global)); sol = solve(prob, Tsit5(),saveat=times_obs,reltol=1e-6, abstol=1e-6,save_everystep=false); Array(sol)[2, :]")
    
    return(vals)
  }
  
  negloglik <- function(log_alpha, log_vmax, theta_c,log_d, n0, log_sigma){
    
    alpha <- exp(log_alpha)
    vmax  <- exp(log_vmax)
    cpar  <- 1 / (1 + exp(-theta_c))
    dpar  <- exp(log_d)
    sigma <- exp(log_sigma)
    r0 <- 10 # fixed arbitrarily
    
    nvals <- satdecay.ode.local(x, alpha, vmax, cpar, dpar, r0, n0)
    
    if (any(!is.finite(nvals))) return(Inf)
    
    -sum(dnorm(y, mean = nvals, sd = sigma, log = TRUE))
  }
  
  # try first fit
  fit.satdecay.ode <- try(mle2(
    negloglik,
    start = list(
      log_alpha = log(alpha.guess),
      log_vmax  = log(vmax.guess),
      theta_c   = qlogis(0.2),
      log_d     = log(d.guess),
      n0   = n0.guess,
      log_sigma = log(0.1)
    ),
    method = "Nelder-Mead",
    control=list(maxit=10000),
    data=data
  ),silent=TRUE)
  if(class(fit.satdecay.ode)=='try-error'){
    print("first attempt at fit failed in get.gr.satdecay.ode") 
    #try again with different settings?
  }
  if(class(fit.satdecay.ode)=='try-error'){ # failed again
    if(!grepl(attr(fit.satdecay.ode,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.satdecay.ode,"condition"))
    }
    #print('fit.satdecay.ode failed after two tries')
  }else{ # can generate plot
    
    # back transform coefficients
    tcoef<-function(cfs){
      vec<-c(exp(cfs[1]),exp(cfs[2]),1 / (1 + exp(-cfs[3])),exp(cfs[4]),cfs[5],exp(cfs[6]))
      names(vec)<-c('alpha','vmax','c','d','n0','sigma')
      vec
    }
    cfs<-data.frame(t(tcoef(coef(fit.satdecay.ode))))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay.ode(x,cfs$alpha,cfs$vmax,cfs$c,cfs$d,10,cfs$n0),min(x),max(x),add=TRUE,col='blue')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay.ode(x,cfs$alpha,cfs$vmax,cfs$c,cfs$d,10,cfs$n0),min(x),max(x),add=TRUE,col='blue')
      }
    }
  }
  
  return(fit.satdecay.ode)
}


#' Extract exponential growth rate assuming exponential death that hits a floor
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances decrease linearly at first, but then hit a floor
#' and remain constant - consistent with a population that declines to the detection
#' limit of fluoresence.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @import zoo
get.gr.flr<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  slopes <- zoo::rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  a.guess<-coef(stats::lm(y~x))[[1]]
  
  fit.flr<-try(nlsLM(y ~ flr(x,a,b,B2,s=1E-10),
                 start=c(a=a.guess,b=min(c(-0.1,round(min(slopes),5))),B2=mean(x)),
                 data = data,
                 upper = c(a=Inf,b=-0.00000001,B2=Inf),
                 control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  
  if(class(fit.flr)=='try-error'){
    if(!grepl(attr(fit.flr,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.flr,"condition"))
    }
    #print('fit.flr failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.flr)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),cfs$B2,max(x),add=TRUE,col='blue')
        graphics::curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),cfs$B2,max(x),add=TRUE,col='blue')
        graphics::curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=TRUE,col='red')
      }
    }
  }
  return(fit.flr)
}


#' Extract exponential growth rate assuming lagged exponential growth that saturates
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances are nearly constant for several time points, before 
#' exponential growth kicks in; subsequently, growth saturates and abundances become 
#' constant again.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @export
#' @importFrom minpack.lm nlsLM nls.lm.control
get.gr.lagsat<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  #slopes <- zoo::rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  #log(max(slopes))
  #a.guess<-coef(lm(y~x))[[1]]
  # round(log(max(slopes)),5)
  
  fit.lagsat<-try(nlsLM(y ~ lagsat(x,a,b,B1,B2,s=1E-10),
                    start = c(B1=mean(x)-(mean(x)-min(x))/2,B2=mean(x)+(max(x)-mean(x))/2, a=min(y)+0.1, b=1),
                    data = data,
                    lower = c(B1=-Inf,B2=-Inf,a=-Inf,b=0.0001),
                    control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  if(class(fit.lagsat)=='try-error'){
    if(!grepl(attr(fit.lagsat,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.lagsat,"condition"))
    }
    #print('fit.lagsat failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.lagsat)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),cfs$B1,cfs$B2,add=TRUE,col='red')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=TRUE,col='blue')
        graphics::curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),cfs$B1,cfs$B2,add=TRUE,col='red')
      }
    }
  }
  return(fit.lagsat)
}


#' Local Slope function
#' 
#' Helper function to calculate and extract the slope of a basic linear regression 
#' relating y to x; the resulting value is used to obtain a reasonable starting guess
#' for the slopes of the piecewise linear functions in \code{lag}, \code{sat}, and
#'  \code{lagsat}
#' 
#' @param d A data frame containing two columns, x and y
#' 
#' @return Slope of the linear regression
#' 
#' @export
localslope<-function (d) {
  m <- stats::lm(y~x, as.data.frame(d))
  return(coef(m)[2])
}

# These two functions help to trim data sets to exclude temperature treatments that are well beyond a species' thermal niche, by identifying the minimum and maximum temperatures where positive growth was observed, and determining the next highest (lowest) temperature treatment.
#funky.low<-function(x,y){
#  dat<-data.frame(x,y) 
#  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
#  Ts<-sort(unique(dat2$x))
  
#  minT<-min(dat2$x[dat2$max.y>0])
#  res<-Ts[max(c(1,which(Ts==minT)-1))]
  
#  return(res)
#}

#funky.high<-function(x,y){
#  dat<-data.frame(x,y) 
#  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
#  Ts<-sort(unique(dat2$x))
#  
#  maxT<-max(dat2$x[dat2$max.y>0])
#  res<-Ts[min(c(length(Ts),which(Ts==maxT)+1))]
#  
#  return(res)
#}
