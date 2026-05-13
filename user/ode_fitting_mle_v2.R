

library(ggplot2)
library(dplyr)
library(growthTools)
library(deSolve)
library(bbmle)
library(tidyr)

expit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}	


## Need some test data

dat<-read.csv("./user/Chan_time_series.csv")
head(dat)

# subset
tmp<-dat[dat$temperature==12 & dat$bacteria==0 & dat$B12==0,]
tmp$ln.fluor<-log(tmp$value)
head(tmp)


#######

# Develop ODE approach (think about adding to package):

tmp %>% filter(dtime>0.1,rep==5) %>%
  ggplot(aes(x=dtime,y=log(value)))+
  geom_point()+
  facet_wrap(~rep)

tmp2<-tmp %>% filter(dtime>0.1,rep==5)
tmp2<-tmp2[order(tmp2$dtime),]

# Differential equation set up
oneRmodel<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dR<- c*q*m*N - (vmax/q)*(R/(R+K))*N*q
    dN<- (vmax/q)*(R/(R+K))*N - m*N
    
    # return the rate of change
    list(c(dR,dN))
  })	# end of with(as.list...
}

times<-seq(0,50,0.01)
N0<-0.01
R0<-10
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)
state<-c(R=R0,N=N0)


# Solve ODE
out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

# Look at results
head(out)

plot(log(N)~time,data=out,type='l')
plot(R~time,data=out,type='l',ylim=c(0,1.05*R0))


# ballpark it:

times<-seq(0,50,0.01)
N0<-exp(2)
R0<-exp(4.5)
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)
state<-c(R=R0,N=N0)

# Solve ODE
out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

plot(log(value)~dtime,data=tmp2,ylim=c(-5,4))
lines(log(N)~time,data=out)



#### Try setting up some fitting:
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)

expit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}	

oneRmodel.ode.NLL<-function(N0,R0,c,q,m,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),q=q,m=m,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  
  # differential equation set up
  oneRmodel<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- c*q*m*N - (vmax/q)*(R/(R+K))*N*q
      dN<- (vmax/q)*(R/(R+K))*N - m*N
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)
  
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  #plot(pops~times,data=pop.data)
  #lines(out[,2]~out[,1])
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(log(tmp2$value),mean=log(preds[,3]),sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  return(NLL)
}

# returns numerical results...
oneRmodel.ode.NLL(N0=exp(2),R0=log(10),c=0.5,q=2,m=0.2,vmax=1,K=1,s=1)

# ballpark initial values:
N0<-exp(2)
R0<-exp(4.5)
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)
state<-c(R=R0,N=N0)

parameters<-c(c=-1.5,q=0.32255274,m=0.04684592,vmax=0.14,K=0.60775770)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=exp(2),R0=log(15),c=-1.5,q=0.3,m=0.04,vmax=0.14,K=0.6,s=1),method="L-BFGS-B",
                        lower=c(N0=0.00001,R0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)
coef(fit.oneRmodel.ode)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=7.463821,R0=4.259327,c=-1.259629,q=0.085585,m=0.04,vmax=0.152865 ,K=0.047065,s=1),method="L-BFGS-B",
                        lower=c(N0=0.00001,R0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)

#pf<-profile(fit.oneRmodel.ode)


times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
# take values from ode solution as predicted values for the given observation times
preds<-out[out[,1]%in%tmp2$dtime,]

# ballpark it:
times<-seq(0,50,0.01)
N0<-7.22884135
R0<-15
#parameters<-c(c=-2.04092566,q=0.32255274,m=0.04684592,vmax=0.14,K=0.60775770)
#parameters<-c(c=-1.5,q=0.32255274,m=0.04684592,vmax=0.14,K=0.60775770)
#parameters<-coef(fit.oneRmodel.ode)[2:6]
parameters<-coef(fit.oneRmodel.ode)[3:7]
#state<-c(R=R0,N=N0)
state<-coef(fit.oneRmodel.ode)[c(2,1)]
names(state)<-c('R','N')

# Solve ODE
out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

plot(log(value)~dtime,data=tmp2,ylim=c(-5,4))
lines(log(N)~time,data=out)



####### not on log scale:


expit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}	

oneRmodel.ode.NLL<-function(R0,N0,c,q,m,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),q=q,m=m,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  
  # differential equation set up
  oneRmodel<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- c*q*m*N - (vmax/q)*(R/(R+K))*N*q
      dN<- (vmax/q)*(R/(R+K))*N - m*N
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)
  
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  plot(tmp2$value~dtime,data=tmp2)
  lines(N~time,data=preds)
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(tmp2$value,mean=preds[,3],sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  return(NLL)
}

# returns numerical results...
oneRmodel.ode.NLL(R0=log(10),N0=exp(2),c=0.5,q=2,m=0.2,vmax=1,K=1,s=1)

# ballpark initial values:
N0<-exp(2)
R0<-exp(4.5)
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)
state<-c(R=R0,N=N0)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(R0=log(15),N0=exp(2),c=-1.5,q=0.3,m=0.04,vmax=0.14,K=0.6,s=1),method="L-BFGS-B",
                        lower=c(R0=0.00001,N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)
coef(fit.oneRmodel.ode)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(R0=15,N0=7,c=-1.5,q=0.05,m=0.15,vmax=0.03,K=0.52,s=10),method="L-BFGS-B",
                        lower=c(R0=0.00001,N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)
coef(fit.oneRmodel.ode)

# suspect overparameterized



fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(R0=15,N0=7,c=-1.5,q=0.05,m=0.15,vmax=0.03,K=0.52,s=10),
                        data=tmp2)
summary(fit.oneRmodel.ode)
coef(fit.oneRmodel.ode)

times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
# take values from ode solution as predicted values for the given observation times
preds<-out[out[,1]%in%tmp2$dtime,]

# ballpark it:
times<-seq(0,50,0.01)
parameters<-coef(fit.oneRmodel.ode)[3:7]
state<-coef(fit.oneRmodel.ode)[c(2,1)]
names(state)<-c('R','N')

# Solve ODE
out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

plot((value)~dtime,data=tmp2)
lines((N)~time,data=out)


parameters<-c(c=-1.5,q=0.32255274,m=0.04684592,vmax=0.14,K=0.60775770)
state<-c(R=R0,N=N0)

out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

plot((value)~dtime,data=tmp2)
lines((N)~time,data=out)



#########


oneRmodel.ode.NLL<-function(N0,c,q,m,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),q=q,m=m,vmax=vmax,K=K)
  state<-c(R=15,N=N0)
  
  # differential equation set up
  oneRmodel<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- c*q*m*N - (vmax/q)*(R/(R+K))*N*q
      dN<- (vmax/q)*(R/(R+K))*N - m*N
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)
  
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  #plot(tmp2$value~dtime,data=tmp2)
  #lines(N~time,data=preds)
  print(c(parameters,state))
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(tmp2$value,mean=preds[,3],sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  print(NLL)
  return(NLL)
}

# returns numerical results...
oneRmodel.ode.NLL(N0=exp(2),c=0.5,q=2,m=0.2,vmax=1,K=1,s=1)

# ballpark initial values:
N0<-exp(2)
R0<-exp(4.5)
parameters<-c(c=0.5,q=2,m=0.2,vmax=1,K=1)
state<-c(R=R0,N=N0)

library(bbmle)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=7,c=-1.5,q=0.05,m=0.15,vmax=0.03,K=0.52,s=10),method="L-BFGS-B",
                        lower=c(N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)
coef(fit.oneRmodel.ode)

# suspect overparameterized

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=7.6041989,c=-1.5485876,q=0.2950981,m=0.1494897,vmax=0.1667206,K=1.5143388,s=0.7355765),method="L-BFGS-B",
                        lower=c(N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
coef(fit.oneRmodel.ode)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=7.0968988,c=-1.5218711,q=0.2934070,m=0.1497804,vmax=0.1749265,K=1.7286851,s=0.7366366),method="L-BFGS-B",
                        lower=c(N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.ode)

fit.oneRmodel.ode<-mle2(oneRmodel.ode.NLL,
                        start=list(N0=7.1563901,c=-1.4782213,q=0.2910410,m=0.1526623,vmax=0.1783534,K=2.0303537,s=0.7133007),method="L-BFGS-B",
                        lower=c(N0=0.00001,c=-100,q=0.00001,m=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)

# inspect the fit:

oneRmodel<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dR<- c*q*m*N - (vmax/q)*(R/(R+K))*N*q
    dN<- (vmax/q)*(R/(R+K))*N - m*N
    
    # return the rate of change
    list(c(dR,dN))
  })	# end of with(as.list...
}

cfs<-coef(fit.oneRmodel.ode)
cfs[2]<-expit(cfs[2])
parameters<-cfs[2:6]
state<-c(R=15,N=cfs[[1]])

out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)

plot((value)~dtime,data=tmp2)
lines((N)~time,data=out)

# looks pretty ok to me...


# how do I extract the maximum?

# helper function to produce N from the ODE system at fixed time ftime:
get.ode.val<-function(ftime){
  times<-unique(sort(c(seq(0,12,0.01),ftime)))
  out<-ode(y=state,times=times,func=oneRmodel,parms=parameters)
  preds<-out[out[,1] %in% ftime,]
  
  return(preds[[3]])
}
#get.ode.val(4)

# use this to run optimization:
opt<-optimize(get.ode.val,c(0,12),maximum=T)
opt

# time where this occurs
opt$maximum

# peak fluorescence
opt$objective

#########

## re-parameterized?


oneRmodel.v2.ode.NLL<-function(R0,N0,c,d,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  
  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)
  
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  #plot(tmp2$value~dtime,data=tmp2)
  #lines(N~time,data=preds)
  print(c(parameters,state))
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(tmp2$value,mean=preds[,3],sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  print(NLL)
  return(NLL)
}

# returns numerical results...
oneRmodel.v2.ode.NLL(R0=15,N0=exp(2),c=0.5,d=0.2,vmax=1,K=1,s=1)


parameters<-c(c=-0.5,d=0.2,vmax=0.5,K=5)
state<-c(R=60,N=7)


# try fitting
fit.oneRmodel.v2.ode<-mle2(oneRmodel.v2.ode.NLL,
                           start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                           lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.v2.ode)
coef(fit.oneRmodel.v2.ode)

fit.oneRmodel.v2.ode<-mle2(oneRmodel.v2.ode.NLL,
                           start=list(R0=58.5086438,N0=6.5805555,c=1.6024987,d=0.5833723,vmax=1.3143025,K=6.0773651,s=1),method="L-BFGS-B",
                           lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),data=tmp2)
summary(fit.oneRmodel.v2.ode)
coef(fit.oneRmodel.v2.ode)



# inspect the fit:

# differential equation set up
oneRmodel.v2<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dR<- vmax*N*(c*d - (R/(R+K)))
    dN<- vmax*N*(R/(R+K)-d)
    
    # return the rate of change
    list(c(dR,dN))
  })	# end of with(as.list...
}


# guesses:
#start=list(R0=15,N0=7,c=-1.5,d=2,vmax=0.03,K=0.52,s=10),method="L-BFGS-B",

parameters<-c(c=-0.5,d=0.2,vmax=0.5,K=5)
state<-c(R=60,N=7)
out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)

plot((value)~dtime,data=tmp2)
lines((N)~time,data=out)


# estimated values
cfs<-coef(fit.oneRmodel.v2.ode)
cfs[3]<-expit(cfs[3])
parameters<-cfs[3:6]
state<-c(R=cfs[[1]],N=cfs[[2]])

out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)

plot((value)~dtime,data=tmp2)
lines((N)~time,data=out)

# looks pretty ok to me...


# how do I extract the maximum?

# helper function to produce N from the ODE system at fixed time ftime:
get.ode.val.v2<-function(ftime){
  times<-unique(sort(c(seq(0,12,0.01),ftime)))
  out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)
  preds<-out[out[,1] %in% ftime,]
  
  return(preds[[3]])
}
#get.ode.val(4)

# use this to run optimization:
opt.v2<-optimize(get.ode.val.v2,c(0,12),maximum=T)
opt.v2

# time where this occurs
opt.v2$maximum

# peak fluorescence
opt.v2$objective

# some improvements from second parameterization
AICtab(fit.oneRmodel.ode,fit.oneRmodel.v2.ode)



#### Get it to run without repeated attempts?

# yes, just bump up the max iterations - not bad:
fit.oneRmodel.v2.ode<-mle2(oneRmodel.v2.ode.NLL,
                           start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                           lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                           control=list(maxit=10000),data=tmp2)
summary(fit.oneRmodel.v2.ode)
coef(fit.oneRmodel.v2.ode)


### Start streamlining this/incorporating into the package?

plot(log(value)~dtime,data=tmp2)
lines(log(N)~time,data=out)

##########


#### Frame this for incorporation into the package ####

# what will I need?

# - needs to work on not-transformed data
# - starting guesses?
# - make sure it's stable
# - needs a get.gr.satdecay() analog: get.gr.satdecay.ode()?
# - integrate into get.gr.rate() functionality
#     - simplify creation of empty data structures? (add new class/object)
# - 




oneRmodel.v2.ode.NLL<-function(R0,N0,c,d,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  
  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)
  
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  #plot(tmp2$value~dtime,data=tmp2)
  #lines(N~time,data=preds)
  #print(c(parameters,state))
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(tmp2$value,mean=preds[,3],sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  #print(NLL)
  return(NLL)
}


# yes, just bump up the max iterations - not bad:
fit.oneRmodel.v2.ode<-mle2(oneRmodel.v2.ode.NLL,
                           start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                           lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                           control=list(maxit=10000),data=tmp2)
summary(fit.oneRmodel.v2.ode)
coef(fit.oneRmodel.v2.ode)

plot(log(value)~dtime,data=tmp2)
lines(log(N)~time,data=out)


### log-transformed version of this:

oneRmodel.log.ode.NLL<-function(R0,N0,c,d,vmax,K,s){
  
  # Generate predicted values
  
  # solve ode given current estimates of parameters
  #R0<-10 # arbitrary...
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  
  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  times<-unique(sort(c(seq(0,12,0.01),tmp2$dtime)))
  
  out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)
  
  # take values from ode solution as predicted values for the given observation times
  preds<-out[out[,1] %in% tmp2$dtime,]
  #print(length(preds))
  
  #plot(tmp2$value~dtime,data=tmp2)
  #lines(N~time,data=preds)
  #print(c(parameters,state))
  
  ### calculate likelihood of these predictions given a normal error distribution
  NLL<- -sum(dnorm(tmp2$ln.fluor,mean=log(preds[,3]),sd=s,log=T))
  
  #	print(c(N0,K,r,NLL))	# when un-commented, this displays function 
  #values being investigated (acts like trace)
  #print(NLL)
  return(NLL)
}


# yes, just bump up the max iterations - not bad:
fit.oneRmodel.log.ode<-mle2(oneRmodel.log.ode.NLL,
                           start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                           lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                           control=list(maxit=10000),data=tmp2)
summary(fit.oneRmodel.log.ode)
coef(fit.oneRmodel.log.ode)

# estimated values
cfs<-coef(fit.oneRmodel.log.ode)
cfs[3]<-expit(cfs[3])
parameters<-cfs[3:6]
state<-c(R=cfs[[1]],N=cfs[[2]])

out<-ode(y=state,times=times,func=oneRmodel.v2,parms=parameters)

plot(ln.fluor~dtime,data=tmp2)
lines(log(N)~time,data=out)



###### Work on packaging up the functions via deSolve approach ####


satdecay.ode<-function(x,R0,N0,c,d,vmax,K,log=F){
  
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  times<-unique(sort(c(seq(0,max(x),0.01),x)))

  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  
  # numerically solve ode
  out<-data.frame(ode(y=state,time=times,func=oneRmodel.v2,parms=parameters))
  
  # take values from ode solution as predicted values for the given observation times
  if(log){
    preds<-log(out$N[out$time %in% x])
  }else{
    preds<-out$N[out$time %in% x]
  }
  
  return(preds)
}

# examples:
satdecay.ode(x=2,R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47)

satdecay.ode(x=c(2,3,4),R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47)

curve(satdecay.ode(x,R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47,log=T),1,11,col='blue',add=T)

# can I avoid the NLL now?

# yes, just bump up the max iterations - not bad:
fit.oneRmodel.log.ode<-mle2(ln.fluor~dnorm(mean=satdecay.ode(tmp2$dtime,R0,N0,c,d,vmax,K,log=T),sd=s),
                            start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                            lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                            control=list(maxit=10000),data=tmp2)

# seems hopeful
summary(fit.oneRmodel.log.ode)
coef(fit.oneRmodel.log.ode)

plot(ln.fluor~dtime,data=tmp2)
curve(satdecay.ode(x,R0 = 54.6,N0=7.17,c=2.58788897,d=0.79855672,vmax=2.17,K=1.47,log=T),1,11,col='blue',add=T)

# cool, now work this into the wrapper function

# (would also be beneficial to optimize code for faster run time)


get.gr.satdecay.ode<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  #slopes <- zoo::rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  
  #a.guess <- coef(stats::lm(y ~ x))[[1]]
  #a.guess <- mean(y[1:2])
  
  #B2.guess<-mean(x) + (max(x) - mean(x))/2
  #B2.guess<-x[which(y==max(y)[1])]
  
  
  # fit.satdecay<-try(nlsLM(y ~ satdecay(x,a,b,b2,B2,s=1E-10),
  #                         start=c(B2=10,a=a.guess,b=round(max(c(slopes,0.0001)),5),
  #                                 b2=round(min(c(slopes,-0.0001)),5)),
  #                         data = data,
  #                         lower = c(B2=-Inf,a=-Inf,b=0.0001,b2=-Inf),
  #                         control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)
  # 
  
  # fitting:
  fit.satdecay.ode<-try(mle2(y~dnorm(mean=satdecay.ode(x,R0,N0,c,d,vmax,K,log=T),sd=s),
                              start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                              lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                              control=list(maxit=10000),data=data))
  if(class(fit.satdecay.ode)=='try-error'){
    fit.satdecay.ode<-try(mle2(y~dnorm(mean=satdecay.ode(x,R0,N0,c,d,vmax,K,log=T),sd=s),
                                    start=list(R0=40,N0=7,c=0.5,d=0.2,vmax=1,K=5,s=10),method="L-BFGS-B",
                                    lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                                    control=list(maxit=10000),data=data))
  }
  if(class(fit.satdecay.ode)=='try-error'){
    if(!grepl(attr(fit.satdecay.ode,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.satdecay.ode,"condition"))
    }
    #print('fit.satdecay failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.satdecay.ode)))
    
    if(plotQ){
      if(!is.na(fpath)){
        grDevices::pdf(fpath)
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay.ode(x,cfs$R0,cfs$N0,cfs$c,cfs$d,cfs$vmax,cfs$K,log=T),min(x),max(x),add=TRUE,col='purple')
        grDevices::dev.off()
      }else{
        graphics::plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        graphics::curve(satdecay.ode(x,cfs$R0,cfs$N0,cfs$c,cfs$d,cfs$vmax,cfs$K,log=T),min(x),max(x),add=TRUE,col='purple')
      }
    }
  }  
  
  return(fit.satdecay.ode)
}

get.gr.satdecay.ode(x=tmp2$dtime,y=tmp2$ln.fluor,plotQ = T)




###### Work on redefining using diffeqr/Julia for speed gains ####

# https://cran.r-project.org/web/packages/diffeqr/readme/README.html
#install.packages("diffeqr")
library(diffeqr)

# installs Julia and additional required packages
#diffeqr::diffeq_setup()



# fails..
de <- diffeqr::diffeq_setup()

f <- function(u,p,t) {
  return(1.01*u)
}

#Then we give it an initial condition and a time span to solve over:
  
u0 <- 1/2
tspan <- c(0., 1.)

#With those pieces we define the ODEProblem and solve the ODE:
  
prob = de$ODEProblem(f, u0, tspan)
sol = de$solve(prob)

#This gives back a solution object for which sol$t are the time points and sol$u are the values. We can treat the solution as a continuous object in time via

sol$.(0.2)

#and a high order interpolation will compute the value at t=0.2. We can check the solution by plotting it:
  
plot(sol$t,sol$u,"l")



### Installation errors?

# [1] "Installed Julia to /Users/colinkremer/Library/Application Support/org.R-project.R/R/JuliaCall/julia/1.9.4/julia-1.9.4"
# Julia version 1.9.4 at location /Users/colinkremer/Library/Application Support/org.R-project.R/R/JuliaCall/julia/1.9.4/julia-1.9.4/bin will be used.
# Loading setup script for JuliaCall...
# LoadError("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/JuliaCall/julia/setup.jl", 16, LoadError("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/JuliaCall/julia/display/RmdJulia.jl", 6, ErrorException("could not load symbol \"SET_SYMVALUE\":\ndlsym(0x91c0ef90, SET_SYMVALUE): symbol not found"))) 
# Error in .julia$cmd(paste0(Rhomeset, "Base.include(Main,\"", system.file("julia/setup.jl",  : 
#                                                                            Error happens when you try to execute command ENV["R_HOME"] = "/Library/Frameworks/R.framework/Resources";Base.include(Main,"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/JuliaCall/julia/setup.jl") in Julia.
#                                                                          To have more helpful error messages,
#                                                                          you could considering running the command in Julia directly
#                                                                          In addition: Warning message:
#                                                                            In system2("bash", "-l -c 'which julia'", stdout = TRUE) :
#                                                                            running command ''bash' -l -c 'which julia'' had status 1


# installed julia directly, via: https://julialang.org/downloads/ and Terminal
# see also https://github.com/JuliaLang/juliaup/wiki/Permission-problems-during-setup
# I found it easier to simply restart Terminal after install than to update path using provided syntax
# BUT - make sure to make note of where Julia gets installed, so you can force julia_setup to use this version


# add older version of Julia
#juliaup add 1.8.3

# try again:
de <- diffeqr::diffeq_setup()


devtools::install_github('SciML/diffeqr', build_vignettes=T)

#install.packages("JuliaCall")
library(JuliaCall)
library(diffeqr)

julia_setup(JULIA_HOME = "/Users/colinkremer/.juliaup/bin",force=T)
#julia_setup(JULIA_HOME = "/Users/colinkremer/.juliaup/bin",version="1.9")

julia_setup(installJulia = T,force=T,version="1.7")

?julia_setup

# try again:
de <- diffeqr::diffeq_setup()

# made it further, but hit new errors:

#' Error: Error happens in Julia.
#' UndefVarError: `DiffEqBase` not defined in `Main`
#' Suggestion: check for spelling errors or missing imports.
#' Hint: DiffEqBase is loaded but not imported in the active module Main.
#' Stacktrace:
#'   [1] getproperty(x::Module, f::Symbol)
#' @ Base ./Base_compiler.jl:47
#' [2] top-level scope
#' @ none:1
#' [3] eval(m::Module, e::Any)
#' @ Core ./boot.jl:489
#' [4] top-level scope
#' @ none:1
#' [5] eval(m::Module, e::Any)
#' @ Core ./boot.jl:489
#' [6] eval_string(x::String)
#' @ Main.JuliaCall /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/JuliaCall/julia/setup.jl:203
#' [7] docall(call1::Ptr{Nothing})
#' @ Main.JuliaCall /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/JuliaCall/julia/setup.jl:0

JuliaCall::autowrap("DiffEqBase.AbstractODESolution", fields = c("t","u"))

?autowrap


#### Seems to be an issue associated with newer versions of R:
# https://github.com/SciML/diffeqr/issues/56

# trying older R (4.4 instead of 4.5) using Rswitch

# reinstalling diffeqr, many dependencies. many many
install.packages("diffeqr")
install.packages("JuliaCall")

library(JuliaCall)
library(diffeqr)

#julia_setup(JULIA_HOME = "/Users/colinkremer/.juliaup/bin",force=T,version="1.10")
julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin",force=T)
#julia_setup(force=T,version="1.10")

# try again:
de <- diffeqr::diffeq_setup()


#####################

# Using R 4.4 arm

# installed Julia 1.10.11 installer from website, determined destination was:
# "/Applications/Julia-1.10.app/Contents/Resources/julia/bin"

# reinstalled JuliaCall, via "install.packages("JuliaCall")"
# install.packages("JuliaCall")

library(JuliaCall)
library(diffeqr)


# then, this worked, after entering Julia via terminal, 
# type: ] to enter package mode
# run: add Pkg
julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin",force=T)

# seems to be running?!?!?
de <- diffeqr::diffeq_setup()

# crap, failed again.

# maybe I need one version older of Julia?



satdecay.ode<-function(x,R0,N0,c,d,vmax,K,log=F){
  
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  times<-unique(sort(c(seq(0,max(x),0.01),x)))
  
  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  
  # numerically solve ode
  out<-data.frame(ode(y=state,time=times,func=oneRmodel.v2,parms=parameters))
  
  # take values from ode solution as predicted values for the given observation times
  if(log){
    preds<-log(out$N[out$time %in% x])
  }else{
    preds<-out$N[out$time %in% x]
  }
  
  return(preds)
}

# examples:
satdecay.ode(x=2,R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47)

satdecay.ode(x=c(2,3,4),R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47)

curve(satdecay.ode(x,R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47,log=T),1,11,col='blue',add=T)

# can I avoid the NLL now?

# yes, just bump up the max iterations - not bad:
fit.oneRmodel.log.ode<-mle2(ln.fluor~dnorm(mean=satdecay.ode(tmp2$dtime,R0,N0,c,d,vmax,K,log=T),sd=s),
                            start=list(R0=60,N0=7,c=-0.5,d=0.2,vmax=0.5,K=5,s=10),method="L-BFGS-B",
                            lower=c(R0=15,N0=0.00001,c=-100,d=0.00001,vmax=0.00001,K=0.00001,s=0.00001),
                            control=list(maxit=10000),data=tmp2)






######### ODE fitting via Julia #########


library(JuliaCall)
library(diffeqr)

# Setup Julia
julia_path <- "/Applications/Julia-1.10.app/Contents/Resources/julia/bin" 
julia_setup(JULIA_HOME = julia_path, installJulia = FALSE)
#julia_command("import DiffEqBase")

de<-diffeqr::diffeq_setup()

# p = (vmax,c,d,k)
f <- function(u,p,t) {
  du1 = p[1]*u[2]*(p[2]*p[3] - u[1]/(u[1]+p[4]))
  du2 = p[1]*u[2]*(u[1]/(u[1]+p[4])-p[3])
  return(c(du1,du2))
}

u0 <- c(54.6,7.17)
tspan <- c(0.0,11.0)
p <- c(2.17,0.93,0.799,1.47)
prob <- de$ODEProblem(f, u0, tspan, p)
sol <- de$solve(prob)

curve(satdecay.ode(x,R0 = 54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47,log=T),1,11,col='blue',add=T)



satdecay.ode<-function(x,R0,N0,c,d,vmax,K,log=F){
  
  parameters<-c(c=expit(c),d=d,vmax=vmax,K=K)
  state<-c(R=R0,N=N0)
  times<-unique(sort(c(seq(0,max(x),0.01),x)))
  
  # differential equation set up
  oneRmodel.v2<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dR<- vmax*N*(c*d - (R/(R+K)))
      dN<- vmax*N*(R/(R+K)-d)
      
      # return the rate of change
      list(c(dR,dN))
    })	# end of with(as.list...
  }
  
  # numerically solve ode
  out<-data.frame(ode(y=state,time=times,func=oneRmodel.v2,parms=parameters))
  
  # take values from ode solution as predicted values for the given observation times
  if(log){
    preds<-log(out$N[out$time %in% x])
  }else{
    preds<-out$N[out$time %in% x]
  }
  
  return(preds)
}


u0 <- c(1.0,0.0,0.0)
tspan <- c(0.0,100.0)
p <- c(10.0,28.0,8/3)
prob <- de$ODEProblem(f, u0, tspan, p)
sol <- de$solve(prob)

mat <- sapply(sol$u,identity)

udf <- as.data.frame(t(mat))

matplot(sol$t,udf,"l",col=1:3)





# p = (vmax,c,d,k)
f <- function(u,p,t) {
  du1 = p[1]*u[2]*(p[2]*p[3] - u[1]/(u[1]+p[4]))
  du2 = p[1]*u[2]*(u[1]/(u[1]+p[4])-p[3])
  return(c(du1,du2))
}

u0 <- c(54.6,7.17)
tspan <- c(0.0,11.0)
p <- c(2.17,0.93,0.799,1.47)
prob <- de$ODEProblem(f, u0, tspan, p)
sol <- de$solve(prob)



julia_command('using DifferentialEquations; f(u,p,t) = -u; u0 = 1.0; tspan = (0.0, 1.0); prob = ODEProblem(f,u0,tspan); sol = solve(prob);')


# Pull the solution back into R
sol <- julia_eval("sol")
print(sol)



julia_eval("using DifferentialEquations")

# define the model in Julia
julia_eval("
function f!(du, u, p, t)
    vmax, c, d, k = p
    du[1] = vmax * u[2] * (c*d - u[1]/(u[1] + k))
    du[2] = vmax * u[2] * (u[1]/(u[1] + k) - d)
end
")

# initial conditions, parameters, time span
julia_eval("u0 = [54.6, 7.17]")
julia_eval("p = [2.17, 0.93, 0.799, 1.47]")
julia_eval("tspan = (0.0, 11.0)")

# define and solve problem
julia_eval("prob = ODEProblem(f!, u0, tspan, p)")
julia_eval("sol = solve(prob, saveat=0:0.1:11.0)")

# bring solution back to R
sol <- julia_eval("Array(sol)")
times <- julia_eval("sol.t") # time vector

df <- data.frame(
  time = times,
  u1 = sol[1, ],
  u2 = sol[2, ]
)

library(ggplot2)

ggplot(df,aes(x=time))+
  geom_line(aes(y=u1),color='blue')+
  geom_line(aes(y=u2),color='green')+
  theme_bw()


## pack this up into a function

satdecay.ode.j<-function(x,R0,N0,c,d,vmax,K,log=F){

  julia_eval("
    function f!(du, u, p, t)
    vmax, c, d, k = p
    du[1] = vmax * u[2] * (c*d - u[1]/(u[1] + k))
    du[2] = vmax * u[2] * (u[1]/(u[1] + k) - d)
    end
    ")
  
  # initial conditions, parameters, time span
  #j.pars<-paste("p = [",vmax,c,d,K,"]",sep="")
  #j.ics<-paste("u0 = [",R0,N0,"]",sep="")
  #j.tspan<-paste("tspan = (",0,max(x),")",sep="")
  
  #print(c(N0,R0,vmax,expit(c),d,K))
  
  julia_assign("u0",c(R0,N0))
  julia_assign("p",c(vmax,c,d,K))
  julia_assign("tspan",c(0,max(x)))
  
  # assign specific times where we want solution values:
  julia_assign("times", x)
  
  #julia_eval("u0 = [54.6, 7.17]")
  #julia_eval("p = [2.17, 0.93, 0.799, 1.47]")
  #julia_eval("tspan = (0.0, 11.0)")
  #julia_eval(j.ics)
  #julia_eval(j.pars)
  #julia_eval(j.tspan)
  
  # define and solve problem
  julia_eval("prob = ODEProblem(f!, u0, tspan, p)")
  julia_eval("sol = solve(prob, saveat=times)")
  
  # bring solution back to R
  sol <- julia_eval("Array(sol)")
  pred_u2 <- sol[2, ]
  
  # guard against solver failure
  if (any(!is.finite(pred_u2))) return(Inf)
  
  if(log){
    pred_u2<-log(pred_u2)
  }
  
  return(pred_u2)
}

# single x (not sure why two values are returned - time = 0 maybe?)
satdecay.ode.j(0.1,R0 = 54.6,N0=7.17,c=expit(0.93),d=0.799,vmax=2.17,K=1.47,log=F)


# multiple x points:
satdecay.ode.j(c(0.1,0.2,0.3),R0 = 54.6,N0=7.17,c=expit(0.93),d=0.799,vmax=2.17,K=1.47,log=F)
satdecay.ode.j(c(0.1,0.2,0.3),R0 = 54.6,N0=7.17,c=expit(0.93),d=0.799,vmax=2.17,K=1.47,log=T)

# plotting
curve(satdecay.ode.j(x,R0 = 54.6,N0=7.17,c=expit(0.93),d=0.799,vmax=2.17,K=1.47,log=F),0.1,11)
curve(satdecay.ode.j(x,R0 = 54.6,N0=7.17,c=expit(0.93),d=0.799,vmax=2.17,K=1.47,log=T),0.1,11)


# run with mle2?
m1<-mle2(ln.fluor~dnorm(mean=satdecay.ode.j(dtime,R0=R0,N0=N0,c=expit(c),d=d,vmax=vmax,K=K,log=T),
                        sd=exp(sigma)),
         start=list(R0=54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47,sigma=0),data=tmp2)


# link functions on pars
m2<-mle2(ln.fluor~dnorm(mean=satdecay.ode.j(dtime,R0=exp(R0),N0=exp(N0),c=expit(c),d=exp(d),vmax=exp(vmax),K=exp(K),log=T),
                        sd=exp(sigma)),
         start=list(R0=log(54.6),N0=log(7.17),c=0,d=log(0.799),vmax=log(2.17),K=log(1.47),sigma=log(1)),data=tmp2)

summary(m2)
# worked reasonably well.

cfs<-data.frame(t(coef(m2)))
cfs$R0

m3<-mle2(ln.fluor~dnorm(mean=satdecay.ode.j(dtime,R0=exp(R0),N0=exp(N0),c=expit(c),d=exp(d),vmax=exp(vmax),K=exp(K),log=T),
                        sd=exp(sigma)),
         start=list(R0=cfs$R0,N0=cfs$N0,c=cfs$c,d=cfs$d,vmax=cfs$vmax,K=cfs$K,sigma=cfs$sigma),data=tmp2)
summary(m3)


# fix R0

m2b<-mle2(ln.fluor~dnorm(mean=satdecay.ode.j(dtime,R0=54.6,N0=exp(N0),c=expit(c),d=exp(d),vmax=exp(vmax),K=exp(K),log=T),
                        sd=exp(sigma)),
         start=list(N0=log(7.17),c=0,d=log(0.799),vmax=log(2.17),K=log(1.47),sigma=log(1)),data=tmp2)

summary(m2b)



# R0 is drifting some. Maybe scale against K? check this is a good fit first, perhaps.

plot(ln.fluor~dtime,data=tmp2)
curve(satdecay.ode.j(x,R0=exp(1.2206e+01),N0=exp(3.2335),c=expit(1.2741e+01),d=exp(5.1406e+00),vmax=exp(-2.0601e+01),K=exp(1.1955),log=T),0.1,11,add=T,col='red')
curve(satdecay.ode.j(x,R0=54.6,N0=7.17,c=0.93,d=0.799,vmax=2.17,K=1.47,log=T),0.1,11,add=T,col='blue')
# no bueno.


###########



library(JuliaCall)
library(bbmle)

julia_setup()
julia_eval("using DifferentialEquations")

# Julia model using alpha
julia_eval("
    function f!(du, u, p, t)
    alpha, vmax, c, d = p
    x = u[1]              # x = u1/k
    u2 = u[2]

    du[1] = alpha * u2 * (c*d - x/(x+1))
    du[2] = vmax * u2 * (x/(x+1) - d)
    end
    ")

#julia_eval("function f!(du,u,p,t); alpha,vmax,c,d = p; x=u[1]; u2=u[2]; du[1]=alpha*u2*(c*d - x/(x+1)); du[2]=vmax*u2*(x/(x+1) - d); end")

# data (replace with yours)
#times_obs <- c(0.0, 0.5, 1.3, 2.7, 5.0, 8.2, 11.0)
#u2_obs <- c(7.17, 6.8, 6.2, 5.5, 4.0, 2.5, 1.8)
times_obs <- tmp2$dtime
u2_obs <- tmp2$ln.fluor

# initial conditions:
# x0 = u1(0)/k → we need k = vmax/alpha
u10 <- 54.6
u20 <- 7.17

tspan <- c(0.0, 11.0)

negloglik <- function(log_alpha, log_vmax, c, log_d, log_sigma) {
  
  # exponentiate (keeps parameters positive)
  alpha <- exp(log_alpha)
  vmax  <- exp(log_vmax)
  c     <- expit(c)
  d     <- exp(log_d)
  sigma <- exp(log_sigma)
  
  # recover k
  k <- vmax / alpha
  
  # initial condition in transformed variable
  x0 <- u10 / k
  
  # send to Julia
  julia_assign("p", c(alpha, vmax, c, d))
  julia_assign("u0", c(x0, u20))
  julia_assign("times", times_obs)
  julia_assign("tspan", tspan)
  
  # solve ODE
  julia_eval('prob = ODEProblem(f!, u0, (tspan[1], tspan[2]), p)')
  julia_eval('sol = solve(prob, Tsit5(), saveat=times,reltol=1e-6, abstol=1e-6)')
  
  sol_mat <- julia_eval("Array(sol)")
  pred_u2 <- sol_mat[2, ]
  
  # guard against solver failure
  if (any(!is.finite(pred_u2))) return(Inf)
  
  # negative log-likelihood
  -sum(dnorm(u2_obs, mean = log(pred_u2), sd = sigma, log = TRUE))
}

# fit
fit <- mle2(
  negloglik,
  start = list(
    log_alpha = log(1.47),
    log_vmax  = log(2.17),
    c     = 0,
    log_d     = log(0.799),
    log_sigma = log(1)
  ),
  method = "Nelder-Mead",
  control=list(maxit=10000)
)

# ~4 seconds


curve(satdecay.ode.j(x,c=0.93,d=,vmax=2.17,K=1.47,log=T),0.1,11,add=T,col='blue')

summary(fit)


# exponentiate (keeps parameters positive)
coef_est <- coef(fit)

alpha <- exp(coef_est["log_alpha"])
vmax  <- exp(coef_est["log_vmax"])
cpar  <- expit(coef_est["c"])
dpar  <- exp(coef_est["log_d"])

# recover k
k <- vmax / alpha

# initial condition in transformed variable
x0 <- u10 / k

# send to Julia
julia_assign("p", c(alpha, vmax, cpar, dpar))
julia_assign("u0", c(x0, u20))
times_plot <- seq(0, 11, by = 0.1)
julia_assign("times", times_plot)
julia_assign("tspan", tspan)

# solve ODE
julia_eval('prob = ODEProblem(f!, u0, (tspan[1], tspan[2]), p)')
julia_eval('sol = solve(prob, Tsit5(), saveat=times,reltol=1e-6, abstol=1e-6)')

sol <- julia_eval("Array(sol)")
t   <- julia_eval("sol.t")

df_fit <- data.frame(
  time = t,
  u2 = log(sol[2, ])
)

plot(df_fit$time, df_fit$u2, type = "l", col = "blue",
     xlab = "Time", ylab = "u2",
     main = "MLE fit")

points(times_obs, u2_obs, col = "red", pch = 16)

# not bad.

# What have I learned?


##########


# Julia model using alpha
julia_eval("
    function f!(du, u, p, t)
    alpha, vmax, c, d = p
    x = u[1]              # x = u1/k
    u2 = u[2]

    du[1] = alpha * u2 * (c*d - x/(x+1))
    du[2] = vmax * u2 * (x/(x+1) - d)
    end
    ")

#julia_eval("function f!(du,u,p,t); alpha,vmax,c,d = p; x=u[1]; u2=u[2]; du[1]=alpha*u2*(c*d - x/(x+1)); du[2]=vmax*u2*(x/(x+1) - d); end")

# data (replace with yours)
#times_obs <- c(0.0, 0.5, 1.3, 2.7, 5.0, 8.2, 11.0)
#u2_obs <- c(7.17, 6.8, 6.2, 5.5, 4.0, 2.5, 1.8)
times_obs <- tmp2$dtime
u2_obs <- tmp2$ln.fluor

# initial conditions:
# x0 = u1(0)/k → we need k = vmax/alpha
u10 <- 54.6 # use fixed value (arbitrary)
#u20 <- 7.17 # to be estimated

tspan <- c(0.0, 11.0)


negloglik <- function(log_alpha, log_vmax, theta_c, log_d, log_sigma, log_u20) {
  
  # transforms
  alpha <- exp(log_alpha)
  vmax  <- exp(log_vmax)
  cpar  <- 1 / (1 + exp(-theta_c))   # logit -> (0,1)
  dpar  <- exp(log_d)
  sigma <- exp(log_sigma)
  u20   <- exp(log_u20)
  
  # recover k and initial x
  k  <- vmax / alpha
  x0 <- u10 / k   # u10 fixed
  
  # send to Julia
  julia_assign("p", c(alpha, vmax, cpar, dpar))
  julia_assign("u0", c(x0, u20))
  julia_assign("times", times_obs)
  
  # solve
  julia_eval("prob = ODEProblem(f!, u0, (0.0,11.0), p)")
  julia_eval("sol = solve(prob, Tsit5(), saveat=times, reltol=1e-6, abstol=1e-6)")
  
  sol_mat <- julia_eval("Array(sol)")
  pred_u2 <- sol_mat[2, ]
  
  if (any(!is.finite(pred_u2))) return(Inf)
  
  -sum(dnorm(u2_obs, mean = log(pred_u2), sd = sigma, log = TRUE))
}

fit2 <- mle2(
  negloglik,
  start = list(
    log_alpha = log(1),
    log_vmax  = log(2),
    theta_c   = qlogis(0.5),  # initial guess c = 0.5
    log_d     = log(0.8),
    log_sigma = log(1),
    log_u20   = log(7)
  ),
  method = "Nelder-Mead",
  control=list(maxit=10000)
)

# didn't need to extend maxit this time?

summary(fit2)

AICtab(fit, fit2)

coef_est <- coef(fit2)

alpha <- exp(coef_est["log_alpha"])
vmax  <- exp(coef_est["log_vmax"])
cpar  <- 1 / (1 + exp(-coef_est["theta_c"]))
dpar  <- exp(coef_est["log_d"])
u20   <- exp(coef_est["log_u20"])

# recover k
k <- vmax / alpha

# initial condition in transformed variable
x0 <- u10 / k

# send to Julia
julia_assign("p", c(alpha, vmax, cpar, dpar))
julia_assign("u0", c(x0, u20))
times_plot <- seq(0, 11, by = 0.1)
julia_assign("times", times_plot)
julia_assign("tspan", tspan)

# solve ODE
julia_eval('prob = ODEProblem(f!, u0, (tspan[1], tspan[2]), p)')
julia_eval('sol = solve(prob, Tsit5(), saveat=times,reltol=1e-6, abstol=1e-6)')

sol2 <- julia_eval("Array(sol)")
t   <- julia_eval("sol.t")

df_fit2 <- data.frame(
  time = t,
  u2 = log(sol2[2, ])
)

plot(ln.fluor~dtime,data=tmp2,xlab = "Time", ylab = "u2",
     main = "MLE fit")
lines(df_fit$time, df_fit$u2, col = "blue")
lines(df_fit2$time, df_fit2$u2, col = "purple")
# so we're seeing run-away value of c here...?


######### making things tidier #########

tspan<-c(0,12)
julia_assign("tspan", tspan)

solve_model_u2 <- function(alpha, vmax, cpar, dpar, u10, u20, times) {
  
  # recover k and initial condition
  k  <- vmax / alpha
  x0 <- u10 / k
  
  # send inputs to Julia
  julia_assign("p", c(alpha, vmax, cpar, dpar))
  julia_assign("u0", c(x0, u20))
  julia_assign("times", times)
  
  # solve (keep this simple and robust)
  julia_eval("prob = ODEProblem(f!, u0, tspan, p)")
  julia_eval("sol = solve(prob, Tsit5(), saveat=times, reltol=1e-6, abstol=1e-6)")
  
  sol_mat <- julia_eval("Array(sol)")
  pred_u2 <- sol_mat[2, ]
  
  return(pred_u2)
}


negloglik <- function(log_alpha, log_vmax, theta_c, log_d, log_sigma, log_u20) {
  
  # transforms
  alpha <- exp(log_alpha)
  vmax  <- exp(log_vmax)
  cpar  <- 1 / (1 + exp(-theta_c))   # logit
  dpar  <- exp(log_d)
  sigma <- exp(log_sigma)
  u20   <- exp(log_u20)
  
  # call black-box solver
  pred_u2 <- solve_model_u2(
    alpha = alpha,
    vmax  = vmax,
    cpar  = cpar,
    dpar  = dpar,
    u10   = u10,
    u20   = u20,
    times = times_obs
  )
  
  # guard
  if (any(!is.finite(pred_u2))) return(Inf)
  
  # likelihood
  -sum(dnorm(log(u2_obs), mean = log(pred_u2), sd = sigma, log = TRUE))
}

negloglik.rawN <- function(log_alpha, log_vmax, theta_c, log_d, log_sigma, log_u20) {
  
  # transforms
  alpha <- exp(log_alpha)
  vmax  <- exp(log_vmax)
  cpar  <- 1 / (1 + exp(-theta_c))   # logit
  dpar  <- exp(log_d)
  sigma <- exp(log_sigma)
  u20   <- exp(log_u20)
  
  # call black-box solver
  pred_u2 <- solve_model_u2(
    alpha = alpha,
    vmax  = vmax,
    cpar  = cpar,
    dpar  = dpar,
    u10   = u10,
    u20   = u20,
    times = times_obs
  )
  
  # guard
  if (any(!is.finite(pred_u2))) return(Inf)
  
  # likelihood
  -sum(dnorm(u2_obs, mean = pred_u2, sd = sigma, log = TRUE))
}


# might try on messier data? or not on log scale


########

times_obs <- tmp2$dtime
u2_obs <- tmp2$value


fit3 <- mle2(
  negloglik,
  start = list(
    log_alpha = log(1),
    log_vmax  = log(2),
    theta_c   = qlogis(0.5),
    log_d     = log(0.8),
    log_sigma = log(1),
    log_u20   = log(7)
  ),
  method = "Nelder-Mead",
  control=list(maxit=10000)
)
summary(fit3)


fit3.raw <- mle2(
  negloglik.rawN,
  start = list(
    log_alpha = log(1),
    log_vmax  = log(2),
    theta_c   = qlogis(0.5),
    log_d     = log(0.8),
    log_sigma = log(1),
    log_u20   = log(7)
  ),
  method = "Nelder-Mead",
  control=list(maxit=10000)
)
summary(fit3.raw)



coef_est <- coef(fit3)

# alpha is getting very large
# defined as vmax/k (suggests that k is getting small)
alpha <- exp(coef_est["log_alpha"])
vmax  <- exp(coef_est["log_vmax"])
cpar  <- 1 / (1 + exp(-coef_est["theta_c"]))
dpar  <- exp(coef_est["log_d"])
u20   <- exp(coef_est["log_u20"])

times_plot <- seq(min(times_obs), max(times_obs), length.out = 300)

pred_u2 <- solve_model_u2(
  alpha = alpha,
  vmax  = vmax,
  cpar  = cpar,
  dpar  = dpar,
  u10   = u10,
  u20   = u20,
  times = times_plot
)

df_plot <- data.frame(
  time = times_plot,
  u2   = log(pred_u2)
)

plot(df_plot$time, df_plot$u2, type = "l", col = "blue",
     lwd = 2,
     xlab = "Time", ylab = "u2",
     main = "MLE fit (black-box solver)")

points(times_obs, log(u2_obs), col = "red", pch = 16)



coef_est <- coef(fit3.raw)

# alpha is getting very large
# defined as vmax/k (suggests that k is getting small)
alpha <- exp(coef_est["log_alpha"])
vmax  <- exp(coef_est["log_vmax"])
cpar  <- 1 / (1 + exp(-coef_est["theta_c"]))
dpar  <- exp(coef_est["log_d"])
u20   <- exp(coef_est["log_u20"])

times_plot <- seq(min(times_obs), max(times_obs), length.out = 300)

pred_u2 <- solve_model_u2(
  alpha = alpha,
  vmax  = vmax,
  cpar  = cpar,
  dpar  = dpar,
  u10   = u10,
  u20   = u20,
  times = times_plot
)

df_plot <- data.frame(
  time = times_plot,
  u2   = pred_u2
)

plot(df_plot$time, df_plot$u2, type = "l", col = "blue",
     lwd = 2,
     xlab = "Time", ylab = "u2",
     main = "MLE fit (black-box solver)")

points(times_obs, u2_obs, col = "red", pch = 16)

# compare estimates of max values

plot(df_plot$time, log(df_plot$u2), type = "l", col = "blue",
     lwd = 2,
     xlab = "Time", ylab = "u2",
     main = "MLE fit (black-box solver)")

points(times_obs, log(u2_obs), col = "red", pch = 16)

# is this really a worse fit than the hard-core near piecewise linear version from before?
# maybe that was just in a tough parameter space.

