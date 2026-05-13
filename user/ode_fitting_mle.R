

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

