source("Some Functions.R")
source("Life cycle function.R")
library(lattice)
library(ggplot2)

Ncycles <- 4
Nyears <- 100
Nboots <- 1000

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

marSurvs <- c(0.077,0.0125)
names(marSurvs) <- c("Good","Bad")

freshSurvMn <- 1
freshSurvCV <- 1e-6
freshRho <- 0

marSurvMn <- 0.075
marSurvCV <- 0.5
marRho <- 0

CR <- 3 # compensation ratio
N0 <- 1e6 # equilibrium adults
propRisk <- 0.485 # fishery closes if population falls below propRisk*N0 

results <- readRDS("results_boot.rds")

par(mar=c(7.5,5,0.5,0.5))
boxplot(UI_risk~age.structure+recCV,data=results,ylab="Closure risk (upper 95% CI)",xlab="",names=NA,col=ifelse(results$recCV=="Historical","dodgerblue","tomato"))
text(1:6, par("usr")[3]-0.1*par("usr")[3],srt = 60, adj= 1, xpd = TRUE,labels = paste(results$age.structure, "Age Str.","\n", results$recCV,"CV",sep=" "), cex=0.75,font=1)
mtext("Age Structure & Variance",side=1,line=5.5,xpd=NA,cex=1)

library(lattice)
bwplot(UI_risk~age.structure|recCV,data=results,outer=TRUE,col="tomato",xlab="Age Structure",ylab="Closure risk (upper 95% CI)",scales=list(cex=1,col="black"),par.settings = list(box.rectangle = list(fill="tomato")))
library(ggplot2)
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
p <- ggplot(data=results, aes(x=age.structure, y=UI_risk, fill=recCV))+ 
  geom_violin(trim=FALSE)+
  facet_wrap(~recCV)+
  geom_jitter(shape=21, bg="grey50",position=position_jitter(0.1),size=2) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Age Structure",y="Closure risk (upper 95% CI)",fill="Variance") +
  theme_minimal() +
  theme(legend.position="top",strip.text.x = element_blank())
p
#ggsave("closure risk.jpeg",units="in",height=5,width=5,dpi=600)
stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")

example <- LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[results$age.structure[37],],freshMarAges=freshMarAges,recCV=recCV[results$recCV[37]],Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk,marSurv_Scen="g-g",probGood=results$probGood[37],probBad=results$probBad[37],goodSurv=marSurvs["Good"],badSurv=marSurvs["Bad"],startProb=0.5)

plot(example$Spawners,example$Recruits,xlab="Spawners",ylab="Recruits")
abline(v=propRisk*N0,lwd=2,col="red")
matplot(example$SpawnersLC,type="l",xlab="Year",ylab="Spawners",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5,ylim=range(c(0,rowSums(example$SpawnersLC))))
lines(rowSums(example$SpawnersLC))
abline(h=propRisk*N0,col="red",lwd=2)
plot(1:Nyears,marSurvs[example$marSurvYrNm],type="l",ylab="Marine survival",xlab="Year")