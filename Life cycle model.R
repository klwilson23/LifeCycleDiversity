source("Some Functions.R")
Ncycles <- 4
Nyears <- 100
CR <- 4 # compensation ratio
N0 <- 1e3/0.2 # equilibrium adults
propRisk <- 0.2 # fishery closes if population falls below propRisk*N0 


freshMarAges <- matrix(NA,ncol=2,nrow=Ncycles)
freshMarAges[1,] <- c(1,2)
freshMarAges[2,] <- c(1,3)
freshMarAges[3,] <- c(2,3)
freshMarAges[4,] <- c(2,2)

lifeCycleNames <- sapply(1:Ncycles,function(x){paste(freshMarAges[x,],collapse=".")})
dimnames(freshMarAges)=list("Life cycles"=lifeCycleNames,
              "Stage"=c("Freshwater","Marine"))

propAge <- rep(1/Ncycles,Ncycles)
propAge <- c(0.4,0.4,0.1,0.1)

recCV <- 1

freshSurvMn <- 0.5
freshSurvCV <- 1e-6
freshRho <- 0.9

marSurvMn <- 0.1
marSurvCV <- 1
marRho <- 0.9

freshSurv <- sapply(1:Ncycles,function(x){exp(-(-log(freshSurvMn)/freshMarAges[x,1]))})

mnfreshSurv <- exp(--log(freshSurvMn)/sum((propAge*freshMarAges[,1])))

marSurv <- sapply(1:Ncycles,function(x){exp(-(-log(marSurvMn)/freshMarAges[x,2]))})

mnmarSurv <- exp(--log(marSurvMn)/sum((propAge*freshMarAges[,2])))

popDyn <- array(NA,dim=c(Nyears,Ncycles,max(freshMarAges[,"Freshwater"])+max(freshMarAges[,"Marine"])),dimnames=list("Year"=1:Nyears,"Life cycles"=lifeCycleNames,"Age"=1:(max(freshMarAges[,"Freshwater"])+max(freshMarAges[,"Marine"]))))

RecruitsLC <- SpawnersLC <- matrix(NA,nrow=Nyears,ncol=Ncycles)
Spawners <- rep(NA,Nyears)
Recruits <- rep(NA,Nyears)
freshSurvYr <- rep(0,Nyears)
marSurvYr <- rep(0,Nyears)

surv.start <- sapply(1:Ncycles,function(x){rep(c(freshSurv[x],marSurv[x]),times=c(freshMarAges[x,1],freshMarAges[x,2]))})
stage.start <- sapply(1:Ncycles,function(x){rep(c(1,2),times=c(freshMarAges[x,1],freshMarAges[x,2]))})
survivorship <- survival <- survivalstage <- matrix(0,nrow=Ncycles,ncol=length(popDyn[1,1,]))
fecundity <- matrix(0,nrow=Ncycles,ncol=length(popDyn[1,1,]))
for(i in 1:Ncycles){
  survival[i,1:length(surv.start[[i]])] <- surv.start[[i]]
  survivalstage[i,1:length(surv.start[[i]])] <- stage.start[[i]]
  survivorship[i,1:length(surv.start[[i]])] <- cumprod(c(1,surv.start[[i]][-length(surv.start[[i]])]))
  fecundity[i,length(surv.start[[i]])] <- 1
}

recProp <- as.vector((t(t(survivorship*fecundity)[t(survivorship*fecundity)!=0]))/sum(survivorship*fecundity))
# in order to get the proportions-at-age observed in spawners, we need recruits to be allocated differentially to those proportions-at-age and survivorship vector
recProp <- (propAge/recProp)/sum((propAge/recProp))
SPR0 <- sum(survivorship*fecundity*recProp) #spawner per recruit
R0 <- N0/SPR0 # equilibrium recruits
alpha.R <- CR/SPR0
beta.R <- -log(alpha.R*SPR0)/(R0*SPR0)

curve(alpha.R*x*exp(beta.R*x),from=0,to=N0,xlab="Spawners",ylab="Recruits")

popDyn[1,,] <- R0*survivorship*recProp
Spawners[1] <- sum(sapply(1:Ncycles,function(x){popDyn[1,x,sum(freshMarAges[x,])]}))
Recruits[1] <- alpha.R*Spawners[1]*exp(beta.R*Spawners[1])

RecruitsLC[1,] <- (alpha.R*Spawners[1]*exp(beta.R*Spawners[1]))*recProp
SpawnersLC[1,] <- sapply(1:Ncycles,function(x){popDyn[1,x,sum(freshMarAges[x,])]})
for(Iyear in 2:Nyears)
{
  # get recruitment from last years' spawners
  Spawners[Iyear] <- max(0,sum(sapply(1:Ncycles,function(x){popDyn[Iyear-1,x,sum(freshMarAges[x,])]})))
  expectedRec <- alpha.R*Spawners[Iyear]*exp(beta.R*Spawners[Iyear])
  Recruits[Iyear] <- rlnorm(1,log(expectedRec),recCV)
  RecruitsLC[Iyear,] <- Recruits[Iyear]*recProp
  SpawnersLC[Iyear,] <- sapply(1:Ncycles,function(x){popDyn[Iyear-1,x,sum(freshMarAges[x,])]})
  
  # calculate time-varying freshwater and marine survival
  freshSurvYr[Iyear] <- rnorm(1,0,freshSurvCV)
  freshSurvYr[Iyear] <- (freshRho)*freshSurvYr[Iyear-1]+(1-freshRho)*freshSurvYr[Iyear]
  #marbetaPars <- get_beta(marSurv,marSurvCV)
  marSurvYr[Iyear] <- rnorm(1,0,marSurvCV)
  marSurvYr[Iyear] <- (marRho)*marSurvYr[Iyear-1]+(1-marRho)*marSurvYr[Iyear]
  for(Ilife in 1:Ncycles)
  {
    popDyn[Iyear,Ilife,1] <- Recruits[Iyear]*recProp[Ilife]
    for(Iage in 2:length(popDyn[Iyear,Ilife,]))
    {
      surv <- ifelse(survivalstage[Ilife,Iage-1]==1,
                     exp(log(freshSurv[Ilife]/(1-freshSurv[Ilife]))+freshSurvYr[Iyear])/(1+exp(log(freshSurv[Ilife]/(1-freshSurv[Ilife]))+freshSurvYr[Iyear])),
                     ifelse(survivalstage[Ilife,Iage-1]==2,exp(log(marSurv[Ilife]/(1-marSurv[Ilife]))+marSurvYr[Iyear])/(1+exp(log(marSurv[Ilife]/(1-marSurv[Ilife]))+marSurvYr[Iyear])),0))
      popDyn[Iyear,Ilife,Iage] <- popDyn[Iyear-1,Ilife,Iage-1]*surv
    }
  }
}
closureRisk <- sum(ifelse(Spawners>(propRisk*N0),0,1))/Nyears

plot(Spawners,Recruits,xlab="Spawners",ylab="Recruits (age 1)")
matplot(SpawnersLC,type="l",xlab="Year",ylab="Spawners",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5)
lines(rowSums(SpawnersLC))

plotMarYr <- exp(log(mnmarSurv/(1-mnmarSurv))+marSurvYr)/(1+exp(log(mnmarSurv/(1-mnmarSurv))+marSurvYr))
plotFreshYr <- exp(log(mnfreshSurv/(1-mnfreshSurv))+freshSurvYr)/(1+exp(log(mnfreshSurv/(1-mnfreshSurv))+freshSurvYr))

plot(plotFreshYr,type="l",ylab="Freshwater survival (yr)",xlab="Year")
plot(plotMarYr,type="l",ylab="Marine survival (yr)",xlab="Year")
