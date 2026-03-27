
# Developing brms fits of growth rate time series

library(dplyr)
library(ggplot2)
library(readxl)

#Raw Counts
counts_PR <- read.csv("/Users/colinkremer/OneDrive - University of Connecticut/Kremer lab/group_research/Mortality/January_2026/Plate Reads/plate_reads_xl/plate_reader_DATA.csv")
#counts_PR <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofConnecticut/Shared Documents - Kremer lab/group_research/Mortality/January_2026/Plate Reads/plate_reads_xl/plate_reader_DATA.csv")
counts_PR <- counts_PR [!(counts_PR$well %in% c("A1", "A2","A3","A4","B1","B2", "B3","B4", "D11", "E11","F11","G11", "H11")),]
counts_PR <-counts_PR [!c(counts_PR $well=="A12"& counts_PR$time==8),]
counts_PR$ID<- paste(counts_PR$well, counts_PR$time)

counts_PR<-counts_PR[!counts_PR$well=="A11",]

average_counts_PR <- counts_PR %>%
  group_by(ID) %>%
  summarize(flour = mean(value, na.rm = TRUE), flour.se = sd(value, na.rm = TRUE)/ sqrt(5))


counts_Flow <- read.csv("/Users/colinkremer/OneDrive - University of Connecticut/Kremer lab/group_research/Mortality/January_2026/Pooled_data.csv")
counts_Flow<-counts_Flow[!c(counts_Flow$well=="A12"& counts_Flow$time==8),]
counts_Flow$wPhyto<-counts_Flow$wLive+counts_Flow$wDead
counts_Flow$well<-paste(counts_Flow$row, counts_Flow$col, sep = "")
counts_Flow$ID<- paste(counts_Flow$well, counts_Flow$time)

df2<- merge(average_counts_PR, counts_Flow, by= "ID" , all.x = TRUE)


PlateID<- read_excel("/Users/colinkremer/OneDrive - University of Connecticut/Kremer lab/group_research/Mortality/January_2026/Plate_Layout.xlsx")
P2<- PlateID[,1:6]
df2$WELL<- df2$well
df<- merge(df2, P2, by= "WELL", all.x = TRUE)
df<- df[df$Acclimated>0,]

df$start.density<- "high"
df$start.density[df$Acute<30]<-"low"

ggplot(data=df, aes(x=wPhyto, y=flour, colour = start.density))+
  geom_errorbar(aes(ymin=flour-flour.se, ymax=flour+flour.se), width=.05, position=position_dodge(0.05))+
  geom_point()+
  #stat_cor(method="pearson", digits = 5)+
  labs(title="wPhyto")+
  theme_bw()


head(df)


ggplot(df,aes(x=time,y=log(wPhyto)))+
  geom_point()+
  geom_line(aes(group=Rep))+
  facet_grid(Acclimated~Acute)

## Pull out a subset:

tmp<-df %>% filter(Acclimated==25,Acute==18)

head(tmp)

ggplot(tmp,aes(x=time,y=log(wPhyto)))+
  geom_point()+
  geom_line(aes(group=Rep))+
  theme_bw()


library(brms)

tmp$ln.N<-log(tmp$wPhyto)
tmp$N<-exp(tmp$ln.N)
tmp$time.days<-tmp$time/24

llogistic<-function(t,K,N0,r){
  log(K/(1+((K-N0)/(N0))*exp(-r*t)))
}

# m.brm<-brm(bf(ln.N ~ llogistic(time,K,N0,r), K~1, N0~1, r~1, nl = TRUE),
#            data=tmp, family = gaussian(),
#            prior=c(prior(normal(8, 1/2), nlpar = "K"),
#                    prior(normal(1, 1/2), nlpar = "N0"),
#                    prior(normal(1, 1/2), nlpar = "r")),
#            control=list(adapt_delta=0.9))


m.brm<-brm(bf(ln.N ~ log(K/(1+((K-N0)/(N0))*exp(-r*time))), K~1, N0~1, r~1, nl = TRUE),
           data=tmp, family = gaussian(),
           prior=c(prior(normal(8, 1/2), nlpar = "K"),
                   prior(normal(1, 1/2), nlpar = "N0"),
                   prior(normal(1, 1/2), nlpar = "r")),
           control=list(adapt_delta=0.9))


plot(wPhyto~time.days,data=tmp)
curve(exp(llogistic(x,2000,50,1.8)),0,5,add=T)

plot(ln.N~time.days,data=tmp)
curve(llogistic(x,2000,50,1.8),0,5,add=T)

m.brm<-brm(bf(ln.N ~ log(K/(1+((K-N0)/(N0))*exp(-r*time))), K~1, N0~1, r~1, nl = TRUE),
           data=tmp, family = gaussian(),
           prior=c(prior(normal(2000, 1/2), nlpar = "K"),
                   prior(normal(50, 1/2), nlpar = "N0"),
                   prior(normal(1.8, 0.2), nlpar = "r")),
           control=list(adapt_delta=0.9))

# incorporate transformations:


# fixing N0 at 50, don't estimate it
# trap K in a link function
m.brm<-brm(bf(ln.N ~ log(exp(Klog)/(1+((exp(Klog)-50)/(50))*exp(-r*time))), Klog~1, r~1, nl = TRUE),
           data=tmp, family = gaussian(),
           prior=c(prior(normal(log(2000), 1/2), nlpar = "Klog"),
                   #prior(normal(50, 1/2), nlpar = "N0"),
                   prior(normal(1.8, 0.2), nlpar = "r")),
           chains=1,
           control=list(adapt_delta=0.9))


# specify initial values for chain:
# https://stackoverflow.com/questions/66050575/specify-initial-values-in-brms
init_func <- function(chain_id=1) {
  list ( Klog = log(2000),
         r  = 1,
         sigma  =  0.8)
}
init_list <- list(
  init_func(chain_id = 1),
  init_func(chain_id = 2)
)


tmp2<-compose_data(tmp)

m.brm<-brm(bf(ln.N ~ log(exp(Klog)/(1+((exp(Klog)-50)/(50))*exp(-r*time))), Klog~1, r~1, nl = TRUE),
           data=tmp2, family = gaussian(),
           prior=c(prior(normal(log(2000), 1/2), nlpar = "Klog"),
                   #prior(normal(50, 1/2), nlpar = "N0"),
                   prior(normal(1.8, 0.2), nlpar = "r")),
           chains=2,
           inits = init_list,
           control=list(adapt_delta=0.9))


# Ok, this at least worked. But the model is not what I want. Was 
# it removing the priors that helped?
m.brm<-brm(bf(ln.N ~ Klog, Klog~1, nl = TRUE),
           data=tmp2, family = gaussian(),
           #prior=c(prior(normal(log(2000), 1/2), nlpar = "Klog")),
           chains=1,
           #inits = init_list,
           control=list(adapt_delta=0.9))

# yes:
m.brm<-brm(bf(ln.N ~ Klog, Klog~1, nl = TRUE),
           data=tmp2, family = gaussian(),
           prior=c(prior(normal(0, 1/2), nlpar = "Klog")),
           chains=1,
           #inits = init_list,
           control=list(adapt_delta=0.9))


# Look at: