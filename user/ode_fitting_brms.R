
### Investigating ODE fitting of abundance time series via brms/Stan ####




#### bayesian version ####

# FOLLOWING: https://www.magesblog.com/post/2021-02-08-fitting-multivariate-ode-models-with-brms/

# we recommend running this in a fresh R session or restarting your current session
#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

library(data.table)
library(brms) 
library(cmdstanr)
library(parallel)
nCores <- detectCores() - 1
options(mc.cores = nCores)

# logit link function:
expit<-function(x){exp(x)/(1+exp(x))}

##### Pull Natalie's data: ########

ndat<-read.csv("./user/Flow_abundance_sample.csv")
head(ndat)

ggplot(ndat,aes(y=wPhyto,x=time))+
  geom_point(aes(color=factor(Sample)))+
  theme_bw()

# formatting
ndat2<-ndat[,c('Sample','time','wPhyto')]
names(ndat2)<-c('Sample','dtime','value')




#### Define model: ###

# ode model in Stan:
OneResource <- "
// Sepcify dynamical system (ODEs)
vector ode_OR(real t, vector y,  vector theta){
  vector[2] dydt;
  
  dydt[1] = theta[1] * y[2] * (theta[2] * theta[3] - (y[1] / (y[1] + theta[4]))); // Resource, R (hare)
  dydt[2] = theta[1] * y[2] * ((y[1] / (y[1] + theta[4])) - theta[3]); // Population, N, (lynx)
  
  return dydt;
}
// Integrate ODEs and prepare output
real OR(real t, real R0, real N0, 
        real vmax, real c, 
        real d, real K){
  vector[2] y0;     // Initial values
  vector[4] theta;  // Parameters
  array[1] vector[2]  y;   // ODE solution
  // Set initial values
  y0[1] = R0; y0[2] = N0;
  // Set parameters
  theta[1] = vmax; theta[2] = c;
  theta[3] = d; theta[4] = K;
  // Solve ODEs
  y = ode_rk45(ode_OR, y0, 0, rep_array(t, 1), theta); 
  // Return relevant population values
  return (y[1,2]);
}
"

# brms syntax for stats model (to be sent to stan)
frmlOR <-  bf(
  value ~ eta,
  nlf(eta ~ log(
    OR(t, R0, N0, vmax, c, d, K)
  )
  ),
  nlf(R0 ~ 2200 * exp(0.3*stdNR0)), # these are transformations to lnormal
  nlf(N0 ~ 90 * exp(0.25*stdNN0)),
  nlf(vmax ~ 2 * exp(0.3 * stdNvmax)),
  #nlf(c ~ 2 * exp(0.3 * stdNc)),
  nlf(c ~ exp(stdNc)/(1+exp(stdNc))), # invokes logit link function, trap 0 < c <1
  nlf(d ~ 0.45 * exp(0.3 * stdNd)),
  nlf(K ~ 4 * exp(0.3 * stdNK)),
  stdNR0 ~ 1,  stdNN0 ~ 1,
  stdNvmax ~ 1, stdNc ~ 1,
  stdNd ~ 1, stdNK ~ 1,
  #sigma ~ 1,
  nl = TRUE)

# define priors (all standard normal distributions)
# This may make the sampling easier, as everything is on a similar scale
# but transformed by the nlf nonlinear transformations above
mypriors <- c(
  prior(normal(0, 1), nlpar = "stdNR0"),
  prior(normal(0, 1), nlpar = "stdNN0"),
  prior(normal(0, 1), nlpar = "stdNvmax"),
  prior(normal(0, 1), nlpar = "stdNc"),
  prior(normal(0, 1), nlpar = "stdNd"),
  prior(normal(0, 1), nlpar = "stdNK")
)

# To build intuition on priors:
#R0 - simulate random numbers from a normal distribution, then transform them:
# hist(2200*exp(0.1*rnorm(1000,0,1)),50)
# hist(2200*exp(0.05*rnorm(10000,0,1)),50)


# more data formatting:

# ode + stan objects to having observations at time = 0
# maybe because of where the initial conditions for the solver are
# implemented. So, shift all time steps by 0.2 to the right (won't 
# change conclusions)
ndat2$t<-round(ndat2$dtime,2)+0.2

# pick out one particular time series to focus on:
ndat3<-ndat2 %>% filter(Sample=="25_B_22")

# run it? hah.
mod <- brm(
  frmlOR, prior = mypriors, 
  stanvars = stanvar(scode = OneResource, block = "functions"),
  data = ndat3, backend = "cmdstan",
  family = brmsfamily("lognormal", link_sigma = "log"),
  control = list(adapt_delta = 0.99),
  seed = 1234, iter = 1000, 
  chains = 4, cores = nCores,
  file = "OneResourceCMDStan_Natalie_25_B_22_trap_c_2.rds")

# Note: if you run more than once, gotta delete the file above or it won't
# overwrite. (Or provide new file name)

# parameter estimates
mod
# (these have to be back-transformed for meaning)

# e.g.,
# pull out vmax chain
rs<-mod %>%
  spread_draws(b_stdNvmax_Intercept)

# 500 left (1000 iterations minus 500 warmup)
range(rs$.iteration)

# on original parameter scale
hist(rs$b_stdNvmax_Intercept,30)

# accounting for back-transformation:
hist(2 * exp(0.3*rs$b_stdNvmax_Intercept),30)


# now do some plotting:

# first time:
expose_functions(mod, vectorize = TRUE, cacheDir = "~/Downloads/")

# data frame for plotting
newdat2<-data.frame(t=seq(0.1,5.1,0.25),dtime=seq(0.1,5.1,0.25),value=NA)

# prediction range for expected trend:
epred <- epred_draws(mod, newdata = newdat2)
head(epred)

ggplot(epred, aes(x = dtime)) +
  stat_lineribbon(aes(y = .epred), 
                  .width = c(.99, .95, .8, .5), 
                  color = "#08519C") +
  #geom_line(data=epred,aes(y='.epred'))+
  #geom_point(data = pred) + 
  #geom_point(data = tmp2) + labs(y="RFU") +
  geom_point(data = ndat3,aes(y=value)) + labs(y="cells") +
  ylim(0,7000)+
  scale_fill_brewer()+
  theme_bw()
# ok, but not mind blowing


# prediction rate for *single observations*
pred <- predicted_draws(mod, newdata = newdat2, ndraws= 1000)

#pred <- predicted_draws(mod, newdata = tmp2, ndraws= 1000)
ggplot(pred, aes(x = dtime, y = value)) +
  stat_lineribbon(aes(y = .prediction), 
                  .width = c(.99, .95, .8, .5), 
                  color = "#08519C") +
  #geom_line(data=epred,aes(y='.epred'))+
  #geom_point(data = pred) + 
  #geom_point(data = tmp2) + labs(y="RFU") +
  geom_point(data = ndat3) + labs(y="cells") +
  ylim(0,15000)+
  scale_fill_brewer()+
  theme_bw()



#### Try with a harder data set:

ggplot(ndat,aes(y=wPhyto,x=time))+
  geom_point(aes(color=factor(Sample)))+
  theme_bw()

# pick out one particular time series to focus on:
ndat3<-ndat2 %>% filter(Sample=="22_A_22")

# run it? hah.
mod <- brm(
  frmlOR, prior = mypriors, 
  stanvars = stanvar(scode = OneResource, block = "functions"),
  data = ndat3, backend = "cmdstan",
  family = brmsfamily("lognormal", link_sigma = "log"),
  control = list(adapt_delta = 0.99),
  seed = 1234, iter = 1000, 
  chains = 4, cores = nCores,
  file = "OneResourceCMDStan_Natalie_22_A_22_trap_c.rds")


#expose_functions(mod, vectorize = TRUE, cacheDir = "~/Downloads/")

# data frame for plotting
newdat2<-data.frame(t=seq(0.1,5.1,0.25),dtime=seq(0.1,5.1,0.25),value=NA)

# prediction range for expected trend:
epred <- epred_draws(mod, newdata = newdat2)
#head(epred)

ggplot(epred, aes(x = dtime)) +
  stat_lineribbon(aes(y = .epred), 
                  .width = c(.99, .95, .8, .5), 
                  color = "#08519C") +
  #geom_line(data=epred,aes(y='.epred'))+
  #geom_point(data = pred) + 
  #geom_point(data = tmp2) + labs(y="RFU") +
  geom_point(data = ndat3,aes(y=value,x=dtime)) + labs(y="cells") +
  ylim(0,7000)+
  scale_fill_brewer()+
  theme_bw()
# ok, but not mind blowing

# pull out vmax chain
rs<-mod %>%
  spread_draws(b_stdNvmax_Intercept)

# 500 left (1000 iterations minus 500 warmup)
range(rs$.iteration)

# on original parameter scale
hist(rs$b_stdNvmax_Intercept,30)

# accounting for back-transformation:
hist(2 * exp(0.3*rs$b_stdNvmax_Intercept),30)

# point estimate for vmax:
median(2 * exp(0.3*rs$b_stdNvmax_Intercept))


