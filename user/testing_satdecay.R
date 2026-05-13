
library(ggplot2)
library(dplyr)
library(growthTools)
library(deSolve)
library(bbmle)
#library(reshape2)
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

# compare growth rate estimates
#bob<-dcast(data.frame(res),rep~method,value.var='mu')
bob<-res %>% pivot_wider(id_cols=rep,names_from=method,values_from=mu)
ggplot(bob,aes(x=no.satdecay,y=all))+
  geom_point()+
  geom_abline()

# compare R2
#bob<-dcast(res,rep~method,value.var='r2')
bob<-res %>% pivot_wider(id_cols=rep,names_from=method,values_from=r2)
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


# On dissecting model objects:

gdat$grs[[1]]
names(gdat$grs[[1]])

gdat$grs[[1]]$models
names(gdat$grs[[1]]$models)

summary(gdat$grs[[1]]$models$gr.lagsat)
coef(gdat$grs[[1]]$models$gr.lagsat)

# here's a specific lagsat fit:
gdat$grs[[1]]$models$gr.lagsat

# save pars for specific lagsat fit:
cfs<-data.frame(t(coef(gdat$grs[[1]]$models$gr.lagsat)))

# estimated saturating density for lagsat model:
lagsat(cfs$B2+0.1,cfs$a,cfs$b,cfs$B1,cfs$B2)

#### 

gdat$grs[[1]]$best.model

# hand one of these to an extractor function:
names(gdat$grs[[1]])

### New functions ####


satval.gr.satdecay<-function(model){
  cfs<-data.frame(t(coef(model)))
  return(cfs$a + cfs$b*cfs$B2)
}
#satval.gr.satdecay(gdat2$grs[[1]]$models$gr.sat)
#names(gdat$grs[[1]]$models)

satval.gr.sat<-function(model){
  cfs<-data.frame(t(coef(model)))
  return(sat(cfs$B2+0.1,cfs$a,cfs$b,cfs$B2))
}

satval.gr.lagsat<-function(model){
  cfs<-data.frame(t(coef(model)))
  return(lagsat(cfs$B2+0.1,cfs$a,cfs$b,cfs$B1,cfs$B2))
}

satval.gr.flr<-function(model){
  cfs<-data.frame(t(coef(model)))
  return(flr(cfs$B2+0.1,cfs$a,cfs$b,cfs$B2))
}

get.longterm.abd<-function(model){
  
  switch(model$best.model,
         gr={satval<-NA},
         gr.sat={satval<-satval.gr.sat(model$models$gr.sat)},
         gr.lag={satval<-NA},
         gr.lagsat={satval<-satval.gr.lagsat(model$models$gr.lagsat)},
         gr.flr={satval<-satval.gr.flr(model$models$gr.flr)},
         gr.satdecay={satval<-satval.gr.satdecay(model$models$gr.satdecay)},
         stop(print("unrecognized growth model type in get.longterm.abd!")))      
  
  return(satval)
}

# single values:
get.longterm.abd(gdat$grs[[1]])
get.longterm.abd(gdat$grs[[3]])

# apply to whole list of fits:
sapply(gdat$grs,get.longterm.abd)
gdat$best.model

# here are the corresponding model types:
sapply(gdat$grs,FUN=function(x) x$best.model)




# need to integrate into package...





