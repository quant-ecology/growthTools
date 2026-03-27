
library(ggplot2)


mgomp<-function(t,topt,d,b0,A,umax,L){
  denom <- (topt^(2))^d
  b0 + (A*exp (-exp (umax*exp (1)/A*(L - t) + 1)))/((((L - t)^(2))^d)/denom + 1)  
}


# 
# b0 = 31;
# A = 126.5;
# umax = 95.8;
# L = 1.31;
# d = 3.19;
# topt = 8.19;

curve(mgomp(x,8.19,3.19,31,126.5,95.8,1.31),0,20)


curve(log(mgomp(x,8.19,3.19,31,126.5,95.8,1.31)),0,20)

## Need some test data

dat<-read.csv("./user/Chan_time_series.csv")
head(dat)

plot(value~dtime,data=dat[dat$temperature==12 & dat$bacteria==0,])

plot(log(value)~dtime,data=dat[dat$temperature==12 & dat$bacteria==0 & dat$B12==0,])

tmp<-dat[dat$temperature==12 & dat$bacteria==0 & dat$B12==0,]

tmp



ggplot(tmp,aes(x=dtime,y=(value)))+
  geom_point()+
  facet_wrap(~rep)

ggplot(tmp,aes(x=dtime,y=log(value)))+
  geom_point()+
  facet_wrap(~rep)

# rep 14 looks pretty good, also 5

tmp1 <- tmp[tmp$rep==14,]
tmp1 <- tmp1[tmp1$dtime>0,]

plot(log(value)~dtime,data=tmp1,xlim=c(0,15),ylim=c(1.5,4))

curve(log(mgomp(x,topt = 5.2,d = 1.3,b0 = exp(1.8),A = 126.5,umax = 15,L = 1.8)),0,20,add=T)

x<-tmp1$dtime
y<-log(tmp1$value)
data<-data.frame(x=x,y=y)

library(minpack.lm)

fit.mgomp<-try(nlsLM(y ~ log(mgomp(x,topt,d,b0,A,umax,L)),
                   start=c(topt = 5.2,d = 1.3,b0 = exp(1.8),A = 126.5,umax = 15,L = 1.8),
                   data = data,
                   #lower = c(B2=-Inf,a=-Inf,b=0.0001),
                   control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)

summary(fit.mgomp)


plot(log(value)~dtime,data=tmp1,xlim=c(0,15),ylim=c(1.5,4))

curve(log(mgomp(x,topt = 5.2,d = 1.3,b0 = exp(1.8),A = 126.5,umax = 15,L = 1.8)),0,20,add=T)
curve(log(mgomp(x,topt = 2.68576,d = 0.99676,b0 = 6.21149,A = 537.63323,umax = 24.64616,L = 4.16433)),0,20,col='red',add=T)

# ok, so the fit fits... is the result interpretable? especially the growth rate... 


fit.mgomp2<-try(nlsLM(exp(y) ~ mgomp(x,topt,d,b0,A,umax,L),
                     start=c(topt = 5.2,d = 1.3,b0 = exp(1.8),A = 126.5,umax = 15,L = 1.8),
                     data = data,
                     #lower = c(B2=-Inf,a=-Inf,b=0.0001),
                     control = nls.lm.control(maxiter=1000,maxfev=1000)),silent=TRUE)

summary(fit.mgomp2)

plot((value)~dtime,data=tmp1)
curve(mgomp(x,topt = 4.0574,d = 0.9228,b0 = 10.6542,A = 91.1532,umax = 20.5685,L = 2.2613),0,20,col='blue',add=T)

# no, this umax growth rate is not interpretable as an exponential growth rate


# try new fitting function:

plot((value)~dtime,data=tmp1)

localslope<-function (d) {
  m <- stats::lm(y~x, as.data.frame(d))
  return(coef(m)[2])
}

get.gr.satdecay(x=tmp1$dtime,y=log(tmp1$value),plotQ =T)




