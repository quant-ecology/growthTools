



library(JuliaCall)
library(bbmle)

julia_setup()
julia_eval("using DifferentialEquations")

# ODE in log-space
julia_eval("
function f!(du,u,p,t)
    alpha, vmax, c, d = p

    x = u[1]
    y = u[2]   # log(u2)

    du[1] = alpha * exp(y) * (c*d - x/(x+1))
    du[2] = vmax * (x/(x+1) - d)
end
")

y_obs <- log(u2_obs)

solve_model_log <- function(alpha, vmax, cpar, dpar, u10, u20, times) {
  
  k  <- vmax / alpha
  x0 <- u10 / k
  y0 <- log(u20)
  
  julia_assign("p", c(alpha, vmax, cpar, dpar))
  julia_assign("u0", c(x0, y0))
  julia_assign("times", times)
  
  julia_eval("prob = ODEProblem(f!, u0, (0.0,11.0), p)")
  julia_eval("sol = solve(prob, Tsit5(), saveat=times,
                          reltol=1e-6, abstol=1e-6)")
  
  sol <- julia_eval("Array(sol)")
  yhat <- sol[2, ]
  
  return(yhat)
}

# another issue to remember: how do I get exp growth rates out of this?

# - see notes on paper, ~= Vmax-d



# things that I want to have:

# satdecay.ode = black-boxed ODE calculation
# get.gr.satdecay.ode = function that receives abundance data, uses mle to fit satdecay.ode
#   - need means of guessing good starting pars
# integration of above model into get.growth.rate function


# satdecay<-function(x,a,b,b2,B2,s=1E-10){
#   a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s) - (-b2/2)*(x-B2) - sqfunc(B2-x,-b2,s)
# }


julia_setup()
julia_eval("using DifferentialEquations")

# ODE in log-space
julia_eval("
function f!(du,u,p,t)
    alpha, vmax, c, d = p

    r = u[1]
    n = u[2]   # log(N)

    du[1] = alpha * exp(n) * (c*d - r/(r+1))
    du[2] = vmax * (r/(r+1) - d)
end
")

# evaluate a dummy/template version of the problem:
julia_eval("u0 = [1.0, 1.0]; p  = [1.0, 1.0, 0.5, 0.5]; tspan = (0.0, 10.0); prob_template = ODEProblem(f!, u0, tspan, p)")

satdecay.ode <- function(x, alpha, vmax, cpar, dpar, r0, n0) {
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

# examples
# multiple values (works)
# satdecay.ode(c(1,2,3),3,1,0.1,0.05,10,1)

# single value (probably busted)
# satdecay.ode(c(1),3,1,0.1,0.05,10,1) # busted:

# need something like this to adapt to single values (rare use case)
#sol <- julia_eval("Array(sol)")
#nvals <- sol[2, ]

# interpolate and return specific values of interest:
# julia_assign("tvals", x)
# if(length(x)==1){
#   vals <- julia_eval("sol(tvals)[2]")
# }else{
#   vals <- julia_eval("[sol(t)[2] for t in tvals]")    
# }



# Now back to MLE fitting

y_obs<-tmp2$ln.fluor
times_obs<-tmp2$dtime

negloglik_log <- function(log_alpha, log_vmax, theta_c,
                          log_d, n0, log_sigma) {
  
  alpha <- exp(log_alpha)
  vmax  <- exp(log_vmax)
  cpar  <- 1 / (1 + exp(-theta_c))
  dpar  <- exp(log_d)
  sigma <- exp(log_sigma)
  #n0   <- exp(log_n0)
  r0 <- 10
  
  nvals <- satdecay.ode(times_obs, alpha, vmax, cpar, dpar, r0, n0)
  
  if (any(!is.finite(nvals))) return(Inf)
  
  -sum(dnorm(y_obs, mean = nvals, sd = sigma, log = TRUE))
}

# mle fit

# precompile solver:
julia_eval("solve(prob_template, Tsit5())")

# run fit
fit_log <- mle2(
  negloglik_log,
  start = list(
    log_alpha = log(0.15),
    log_vmax  = log(0.5),
    theta_c   = qlogis(0.2),
    log_d     = log(0.1),
    n0   = 2,
    log_sigma = log(0.2)
  ),
  method = "Nelder-Mead",
  control=list(maxit=10000)
)

# works, but back to being pretty slow...
# and even slower once I tried to fix n0...
summary(fit_log)

# try optimizing starting guesses by hand:
xs<-seq(0,11,0.1)
ys<-satdecay.ode(xs,0.15,.5,0.2,0.1,10,2)
plot(ln.fluor~dtime,data=tmp2); lines(ys~xs)



# back transform coefficients
tcoef<-function(cfs){
  vec<-c(exp(cfs[1]),exp(cfs[2]),1 / (1 + exp(-cfs[3])),exp(cfs[4]),cfs[5],exp(cfs[6]))
  names(vec)<-c('alpha','vmax','c','d','n0','sigma')
  vec
}
tcoef(coef(fit_log))
cfs<-data.frame(t(tcoef(coef(fit_log))))

# how does it look?

xs<-seq(0,11,0.1)
ys<-satdecay.ode(xs,cfs$alpha,cfs$vmax,cfs$c,cfs$d,10,cfs$n0)
#plot(ln.fluor~dtime,data=tmp2); lines(ys~xs)
lines(ys~xs,col='blue')

# result didn't really change. code wasn't that much faster... but this might be survivable as is.



# start packaging this up:
# could potentially see large gains from moving even more of the likelihood calcs into Julia
get.gr.satdecay.ode<-function(x,y,plotQ=F,fpath=NA,id=''){
  data<-data.frame(x=x,y=y)
  
  if(length(unique(x))<4){
    print("error: fewer than four distinct time steps provided to get.gr.satdecay.ode")
    break;
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
  n0.guess <- y[x==min(x)[1]]
  c.guess <- 0.2
  alpha.guess <- 0.1 * vmax.guess # careful with this one; linked to r0 assumption
  
  # precompile solver: (is this necessary/helpful?)
  julia_eval("prob = remake(prob_template); solve(prob, Tsit5(), saveat=times_obs)")
  
  # set up for likelihood calculation:
  julia_assign("times_obs", x)
  julia_assign("tmax_global", max(max(x), 10))
  
  # local version of satdecay.ode(), to optimize run time. Uses fixed time vals
  satdecay.ode.local <- function(x, alpha, vmax, cpar, dpar, r0, n0) {
    # define time range
    tmax <- max(max(x), 10)
    
    julia_assign("p_new", c(alpha, vmax, cpar, dpar))
    julia_assign("u0_new", c(r0, n0))
    
    # below only works if x is more than one value
    vals <- julia_eval("prob = remake(prob_template,u0=u0_new,p=p_new,tspan=(0.0, tmax_global)); sol = solve(prob, Tsit5(),saveat=times_obs,reltol=1e-6, abstol=1e-6,save_everystep=false); Array(sol)[2, :]")
    
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

# 1.234
system.time(fit.ode<-get.gr.satdecay.ode(tmp2$dtime,tmp2$ln.fluor,plotQ = T,id='new result!'))

# busted
#library(growthTools)
library(minpack.lm)
system.time(fit.reg<-get.gr.satdecay(tmp2$dtime,tmp2$ln.fluor,plotQ = T,id='new result - old method!'))

AICtab(fit.ode,fit.reg)


# how does it look?
summary(fit)
vcov(fit)


# what if we try on other data, how robust is it?

dat<-read.csv("./user/Chan_time_series.csv")
head(dat)

# subset
tmp<-dat[dat$temperature==12 & dat$bacteria==0 & dat$B12==0,]
tmp$ln.fluor<-log(tmp$value)
head(tmp)


tmp %>% filter(dtime>0.1) %>%
  ggplot(aes(x=dtime,y=log(value)))+
  geom_point()+
  facet_wrap(~rep)

# pull diff't rep:
tmp3<-tmp %>% filter(dtime>0.1,rep==3)
tmp3<-tmp3[order(tmp3$dtime),]


system.time(fit<-get.gr.satdecay.ode(tmp3$dtime,tmp3$ln.fluor,plotQ = T,id='rep3!'))


bob<-dat[dat$bacteria==0 & dat$B12==0,]

bob %>% filter(dtime>0.1) %>%
  ggplot(aes(x=dtime,y=log(value)))+
  geom_point()+
  facet_grid(temperature~rep)

bob %>% filter(dtime>0.1) %>%
  ggplot(aes(x=dtime,y=log(value)))+
  geom_point(aes(color=as.factor(rep)))+
  facet_wrap(~temperature)




