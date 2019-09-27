source("Some Functions.R")
source("Life cycle function.R")
library(lattice)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(broom)
library(wesanderson)
Ncycles <- 4
Nyears <- 100
Nboots <- 100

freshMarAges <- matrix(NA,ncol=2,nrow=Ncycles)
freshMarAges[1,] <- c(1,2)
freshMarAges[2,] <- c(1,3)
freshMarAges[3,] <- c(2,2)
freshMarAges[4,] <- c(2,3)

lifeCycleNames <- sapply(1:Ncycles,function(x){paste(c(sum(freshMarAges[x,])+1,freshMarAges[x,1]+1),collapse=".")})
dimnames(freshMarAges)=list("Life cycles"=lifeCycleNames,
                            "Stage"=c("Freshwater","Marine"))

propAge1 <- c(0.545,0.308,0.096,0.051)
propAge2 <- c(0.459,0.497,0.021,0.023)
propAge3 <- rep(1/Ncycles,Ncycles)
propAge <- rbind(propAge1,propAge2,propAge3)
row.names(propAge) <- c("Historical","Recent","Even")
recCV <- c(0.486,0.729)
names(recCV) <- c("Historical","Recent")
probGood <- c(0.2,0.4,0.8)
probBad <- c(0.2,0.4,0.8)
startProb <- 0.5

marSurvs <- c(0.077,1.25*0.0125)
names(marSurvs) <- c("Good","Bad")

freshSurvMn <- 1
freshSurvCV <- 1e-6
freshRho <- 0

marSurvMn <- 0.075
marSurvCV <- 0.5
marRho <- 0

CR <- 3 # compensation ratio
N0 <- 1e6 # equilibrium adults
propRisk <- 0.1#0.485 # fishery closes if population falls below propRisk*N0 

results <- readRDS("results_boot.rds")

library(ggplot2)
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
p <- ggplot(data=results, aes(x=age.structure, y=mn_risk, fill=recCV))+ 
  geom_violin(trim=TRUE)+
  facet_wrap(~recCV)+
  geom_jitter(shape=21, bg="grey50",position=position_jitter(0.1),size=2) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Age Structure",y="Closure risk (mean)",fill="Variance") +
  theme_minimal() +
  theme(legend.position="top",strip.text.x = element_blank())
p
ggsave("closure risk.jpeg",units="in",height=5,width=5,dpi=600)
stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")
results$scen <- factor(paste("G-G ",results$probGood," B-B ",results$probBad,sep=""),
                       levels=c("G-G 0.2 B-B 0.2","G-G 0.4 B-B 0.2","G-G 0.8 B-B 0.2",
                                "G-G 0.2 B-B 0.4","G-G 0.4 B-B 0.4","G-G 0.8 B-B 0.4",
                                "G-G 0.2 B-B 0.8","G-G 0.4 B-B 0.8","G-G 0.8 B-B 0.8"),
                       labels=c("Lo.\nLo.","Med.\nLo.","Hi.\nLo.",
                                "Lo.\nMed","Med.\nMed.","Hi.\nMed",
                                "Lo.\nHi.","Med.\nHi.","Hi.\nHi."))
results$history <- factor(paste(results$age.structure," ages, ",results$recCV, " CV",sep=""))
p <- ggplot(data=results, aes(x=scen, y=mn_risk,colour=history))+
  geom_point(data=results,aes(y=mn_risk,colour=history))+
  geom_linerange(data=results,aes(ymin=LI_risk, ymax=UI_risk,colour=history))+
  facet_wrap(~history,nrow=2,dir="v")+
  scale_fill_brewer(palette="Dark2") +
  labs(x="Correlation severity",y="Closure risk (95% CI)",colour="Scenario") +
  theme_minimal() +
  theme(legend.position="top",strip.text.x = element_blank())
pAnnotated <- annotate_figure(p,bottom=text_grob(wrapper("Closure risks for Skeena sockeye under alternate scenarios of age structure diversity, recruitment variation, and correlations in marine survival conditions (top: Good-to-Good correlation, bottom: Bad-to-Bad correlation).",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("closure risk v2.jpeg",plot=pAnnotated,units="in",height=5,width=8,dpi=600)

example <- LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[results$age.structure[13],],freshMarAges=freshMarAges,recCV=recCV[results$recCV[13]],Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk,marSurv_Scen="g-g",probGood=results$probGood[13],probBad=results$probBad[13],goodSurv=marSurvs["Good"],badSurv=marSurvs["Bad"],startProb=0.5)

example2 <- LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[results$age.structure[25],],freshMarAges=freshMarAges,recCV=recCV[results$recCV[25]],Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk,marSurv_Scen="g-g",probGood=results$probGood[25],probBad=results$probBad[25],goodSurv=marSurvs["Good"],badSurv=marSurvs["Bad"],startProb=0.5)

example3 <- LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[results$age.structure[52],],freshMarAges=freshMarAges,recCV=recCV[results$recCV[52]],Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk,marSurv_Scen="g-g",probGood=results$probGood[52],probBad=results$probBad[52],goodSurv=marSurvs["Good"],badSurv=marSurvs["Bad"],startProb=0.5)


plot(example$Spawners,example$Recruits,xlab="Spawners",ylab="Recruits")
plot(example3$Spawners,example3$Recruits,xlab="Spawners",ylab="Recruits")

abline(v=propRisk*N0,lwd=2,col="red")


jpeg("time series.jpeg",height=7,width=8,units="in",res=800)
layout(1)
par(mar=c(5,5,1,1))
matplot(example$SpawnersLC,type="l",xlab="Year",ylab="Spawners",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5,ylim=range(c(0,rowSums(example$SpawnersLC))))
lines(rowSums(example$SpawnersLC))
abline(h=propRisk*N0,col="red",lwd=2)
dev.off()

jpeg("proportions in time.jpeg",height=7,width=8,units="in",res=800)
layout(1)
par(mar=c(5,5,1,1))
matplot(example$SpawnersLC/rowSums(example$SpawnersLC),type="l",xlab="Year",ylab="Proportion of metapopulation",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5,ylim=c(0,1))
dev.off()

jpeg("survival scenarios.jpeg",height=7,width=8,units="in",res=800)
layout(matrix(1:3,nrow=3))
par(mar=c(0,5,5,18))
plot(1:Nyears,marSurvs[example$marSurvYrNm],type="l",ylab="",xlab="",lwd=2,col=adjustcolor("dodgerblue",0.7),xpd=F)
par(mar=c(2.5,5,2.5,18))
plot(1:Nyears,marSurvs[example2$marSurvYrNm],type="l",lwd=2,col=adjustcolor("orange",0.7),ylab="",xlab="",xpd=F,cex.lab=1.5)
mtext("Marine survival (per ocean life stage)",side=2,line=2.5)
par(mar=c(5,5,0,18))
plot(1:Nyears,marSurvs[example3$marSurvYrNm],type="l",lwd=2,col=adjustcolor("tomato",0.7),ylab="",xlab="",xpd=F,cex.lab=1.5)
mtext("Year",side=1,line=2)
legend(105,0.15,c("G-G (0.8) & B-B (0.2)","G-G (0.8) & B-B (0.4)","G-G (0.8) & B-B (0.8)"),bty="n",lwd=2,col=adjustcolor(c("dodgerblue","orange","tomato")),title="Scenario correlations",xpd=NA,cex=1.5)
dev.off()