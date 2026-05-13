
library(ggplot2)
library(dplyr)

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

# subset
tmp<-dat[dat$temperature==12 & dat$bacteria==0 & dat$B12==0,]
tmp$ln.fluor<-log(tmp$value)
head(tmp)

ggplot(tmp,aes(x=dtime,y=(value)))+
  geom_point()+
  facet_wrap(~rep)

tmp %>% filter(dtime>0.1) %>%
ggplot(aes(x=dtime,y=log(value)))+
  geom_point()+
  facet_wrap(~rep)

## Try out new method

# all methods
gdat <- tmp %>% filter(dtime>0.1) %>% group_by(rep) %>% do(grs=get.growth.rate(.$dtime,.$ln.fluor,.$id,plot.best.Q=T,fpath=NA))
res1<-gdat %>% summarise(rep,mu=grs$best.slope,best.model=grs$best.model,r2=grs$best.model.slope.r2)
res1
# looks like satdecay was pretty commonly best model


# everything but satdecay
gdat2 <- tmp %>% filter(dtime>0.1) %>% group_by(rep) %>% do(grs=get.growth.rate(.$dtime,.$ln.fluor,.$id,plot.best.Q=T,fpath=NA,methods = c("linear", "lag", "sat", "flr", "lagsat")))
res2<-gdat2 %>% summarise(rep,mu=grs$best.slope,best.model=grs$best.model,r2=grs$best.model.slope.r2)
res2

# looks like satdecay was pretty commonly best model

# compare results
res1$method<-"all"
res2$method<-"no.satdecay"

res<-rbind(res1,res2)
head(res)

library(reshape2)

# compare growth rate estimates
bob<-dcast(res,rep~method,value.var='mu')
ggplot(bob,aes(x=no.satdecay,y=all))+
  geom_point()+
  geom_abline()

# compare R2
bob<-dcast(res,rep~method,value.var='r2')
ggplot(bob,aes(x=no.satdecay,y=all))+
  geom_point()+
  geom_abline()


# extract extra traits for sat decay model:

get.satdecay.pars<-function(x){
  cfs<-coef(x$models$gr.satdecay)
}

trlist<-lapply(gdat$grs,FUN = get.satdecay.pars)

satdecay.pars<-data.frame(do.call(rbind, trlist))

# b is the increasing growth rate:
satdecay.pars$b

# b2 is the decreasing growth rate:
satdecay.pars$b2

# B2 is the break point (time when direction of growth changes):
satdecay.pars$B2

# peak abundances
satdecay.pars$a + satdecay.pars$b*satdecay.pars$B2



