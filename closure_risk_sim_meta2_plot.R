################################################

# Closure risk simulation & plot Wild + Enhanced separately
# Nov. 28, 2019 - adding bootstrap for metapopulation within each period

################################################
library(tidyverse)
library(Rmisc) #calculate confidence intervals
library(unikn) 
 
#mean abundances for each period
av_abund_wild <- c(1.73e6, 0.88e6, 4.5e5, 0)
av_abund_total <- c(1.73e6, 0.88e6, 4.5e5, 1.71e6)
meta_abund_wild <- c(1.73e6, 0.88e6, 4.5e5, 0)
meta_abund_total <- c(1.73e6, 0.88e6, 4.5e5, 1.71e6)
av_cv_wild <- c(0.977, 1.005, 0.537, 0)
av_cv_total <- c(0.977, 1.005, 0.537, 0.569)
meta_cv_wild <- c(0.483, 0.677, 0.424, 0)
meta_cv_total <- c(0.483, 0.677, 0.424, 0.517)

av_cv_wild_log <- sqrt(log(1+(av_cv_wild*av_abund_wild)^2/(av_abund_wild^2)))
av_cv_total_log <- sqrt(log(1+(av_cv_total*av_abund_total)^2/(av_abund_total^2)))
meta_cv_wild_log <- sqrt(log(1+(meta_cv_wild*meta_abund_wild)^2/(meta_abund_wild^2)))
meta_cv_total_log <- sqrt(log(1+(meta_cv_total*meta_abund_total)^2/(meta_abund_total^2)))
av_abund_wild_log <- log(av_abund_wild)-av_cv_total_log^2/2
av_abund_total_log <- log(av_abund_total)-av_cv_wild_log^2/2
meta_abund_wild_log <- log(meta_abund_wild)-meta_cv_wild_log^2/2
meta_abund_total_log <- log(meta_abund_total)-meta_cv_total_log^2/2


#name the periods
names(meta_abund_wild) <- names(meta_abund_total) <- names(meta_cv_wild) <- names(meta_cv_total) <- names(av_cv_wild) <- names(av_cv_total) <- Period <- c("1","2","3", "4")

#set closure target, Period, and #iterations
fishery_target <- 1.05e6
Nperiods <- length(unique(Period))
Nboots <- 10000

#create 4 blank matrices that each have columns for Nboots & Nperiods
Ntotal <- closure <- matrix(NA, Nboots, Nperiods, 
                                 dimnames = list("Boots" = 1:Nboots, 
                                 "Period" = c("1","2","3", "4")))

av_closure_wild <- av_Nwild <-  matrix(NA, Nboots, Nperiods, 
                                          dimnames = list("Boots" = 1:Nboots, 
                                                          "Period" = c("1","2","3", "4")))
av_closure_tot <- av_Ntot <- matrix(NA, Nboots, Nperiods, 
                                        dimnames = list("Boots" = 1:Nboots, 
                                                        "Period" = c("1","2","3", "4")))
meta_closure_wild <- meta_Nwild <- matrix(NA, Nboots, Nperiods, 
                                         dimnames = list("Boots" = 1:Nboots, 
                                                         "Period" = c("1","2","3", "4")))
meta_closure_tot <- meta_Ntot <- matrix(NA, Nboots, Nperiods, 
                                 dimnames = list("Boots" = 1:Nboots, 
                                                 "Period" = c("1","2","3", "4")))

#create a blank array to receive the data
Nt_pops <- array(NA, dim = c(Nboots, Nperiods), dimnames = list("Boots"=1:Nboots,
                                                                "Period"= c("1","2","3", "4")))
#run a for-loop for each period
for(p in 1:Nperiods)
{
  for(bI in 1:Nboots) #for each bootstrap iteration
  {
    Ntotal[bI,p] <- sum(Nt_pops[bI,p], na.rm=TRUE)
    closure[bI,p] <- ifelse(Ntotal[bI,p] <= fishery_target,1,0)
    # do the simulation for the metapopulation performance within each time period
    # (1) generate random normal draw for metapopulation abundance given mean N and the meta CV
    meta_Ntot[bI,p] <- max(0,rlnorm(1,meanlog=meta_abund_total_log[p],sdlog=meta_cv_total_log[p]))
    # (1b) - repeat this but for wild abundance/CV only
    meta_Nwild[bI,p] <- max(0,rlnorm(1,meanlog=meta_abund_wild_log[p],sdlog=meta_cv_wild_log[p]))
    # calculate whether the wild-only fishery closed
    meta_closure_wild[bI,p] <- ifelse(meta_Nwild[bI,p] <= fishery_target,1,0)
    # calculate whether the total fishery closed
    meta_closure_tot[bI,p] <- ifelse(meta_Ntot[bI,p] <= fishery_target,1,0)
    #(1) generate random normal draw for metapopulation abundance given mean N and the meta CV
    av_Ntot[bI,p] <- max(0,rlnorm(1,meanlog=av_abund_total_log[p],sdlog=av_cv_total_log[p]))
    # (1b) - repeat this but for wild abundance/CV only
    av_Nwild[bI,p] <- max(0,rlnorm(1,meanlog=av_abund_wild_log[p],sdlog=av_cv_wild_log[p]))
    # calculate whether the wild-only fishery closed
    av_closure_wild[bI,p] <- ifelse(av_Nwild[bI,p] <= fishery_target,1,0)
    # calculate whether the total fishery closed
    av_closure_tot[bI,p] <- ifelse(av_Ntot[bI,p] <= fishery_target,1,0)
  }
}

#Quantify Closure results for each component & time-period
#meta CV wild 
metaClosureWild <- as.data.frame(meta_closure_wild) %>%
  gather(period, closure) %>%
  filter(period == "2") %>% #manually change for each time-period (1-4)
  summarize(avg = mean(closure))

#average CV wild
avClosureWild <- as.data.frame(av_closure_wild) %>%
  gather(period, closure) %>%
  filter(period == "2") %>% #manually change for each time-period
  summarize(avg = mean(closure))

#meta CV total (wild + enhancement)
metaClosureTot <- as.data.frame(meta_closure_tot) %>%
  gather(period, closure) %>%
  filter(period == "4") %>%
  summarize(avg = mean(closure))

#average CV total (wild + enhancement)
avClosureTot <- as.data.frame(av_closure_tot) %>%
  gather(period, closure) %>%
  filter(period == "4") %>%
  summarize(avg = mean(closure))

#create data.frames for Total (wild+enhanced) & Wild components for plotting (1 & 2)
#1. Total abundance (average cv)
avTot <- as.data.frame(av_Ntot) %>%
  gather(period, abundance) %>%
  mutate(component = "averaged") 
#total abundance (meta cv)
metaTot <- as.data.frame(meta_Ntot) %>%
  gather(period, abundance) %>%
  mutate(component = "meta")
#combine
ggTotal <- rbind(avTot, metaTot) %>%
  mutate(abundance = abundance/1000) %>% #reduce abundance to thousands 
  mutate(abundance = round(abundance, digits=0)) %>% #round digits to 0
  filter(period == "4") #filter for recent period b/c only these data required for plot

#2. Wild abundance (average cv)
avWild <- as.data.frame(av_Nwild) %>%
  gather(period, abundance) %>%
  mutate(component = "averaged")
#wild abundance (average cv)
metaWild <- as.data.frame(meta_Nwild) %>%
  gather(period, abundance) %>%
  mutate(component = "meta")
#combine
ggWild <- rbind(avWild, metaWild) %>%
  mutate(abundance = abundance/1000) %>% #reduce abundance to thousands 
  mutate(abundance = round(abundance, digits=0)) %>% #round digits to 0
  filter(period < "4")

#plot
seecol(pal_unikn)
col <- c("blue", "#00A9E0")

ggplot(ggWild, aes(x=period, y=abundance, fill=component)) +
  geom_violin(alpha = 0.9, scale = "width") +
  scale_fill_manual(values=col) +
  theme_classic() +
  theme(legend.position="none") +
  labs(x="Fishery period", y="Sockeye abundance (thousands)") +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12)) +
  geom_violin(data = ggTotal, aes(fill=component)) +
  geom_hline(yintercept=1050, size=1, colour="red") +
  scale_x_discrete(labels=c("1" = "1913-1923", "2" = "1933-1947",
                            "3" = "2010-2017", "4" = "2010-2017"))
  ggsave(file = paste0("Fishery_closure_combined_revised.pdf")) 

#add boxplot to violins - "geom_boxplot(aes(fill=component), position = position_dodge(0.9), width=0.03)"

##################################################

