source("Some Functions.R")
Ncycles <- 4
Nyears <- 100

freshMarAges <- matrix(NA,ncol=2,nrow=Ncycles)
freshMarAges[1,] <- c(1,2)
freshMarAges[2,] <- c(1,3)
freshMarAges[3,] <- c(2,3)
freshMarAges[4,] <- c(2,2)

lifeCycleNames <- sapply(1:Ncycles,function(x){paste(freshMarAges[x,],collapse=".")})
dimnames(freshMarAges)=list("Life cycles"=lifeCycleNames,
              "Stage"=c("Freshwater","Marine"))

propAge <- rep(1/Ncycles,Ncycles)
propAge <- c(0.1,0.1,0.1,0.7)

recCV <- 0.35

freshSurvMn <- 0.1
freshSurv <- exp(--log(freshSurvMn)/sum((propAge*freshMarAges[,1])))
freshSurvCV <- 0.5
freshRho <- 0.01

marSurvMn <- 0.1
marSurv <- exp(--log(marSurvMn)/sum((propAge*freshMarAges[,2])))
marSurvCV <- 0.5
marRho <- 0.9

popDyn <- array(NA,dim=c(Nyears,Ncycles,max(freshMarAges[,"Freshwater"])+max(freshMarAges[,"Marine"])),dimnames=list("Year"=1:Nyears,"Life cycles"=lifeCycleNames,"Age"=1:(max(freshMarAges[,"Freshwater"])+max(freshMarAges[,"Marine"]))))

RecruitsLC <- SpawnersLC <- matrix(NA,nrow=Nyears,ncol=Ncycles)

Spawners <- rep(NA,Nyears)
Recruits <- rep(NA,Nyears)
freshSurvYr <- rep(freshSurv,Nyears)
marSurvYr <- rep(marSurv,Nyears)

surv.start <- sapply(1:Ncycles,function(x){rep(c(freshSurv,marSurv),times=c(freshMarAges[x,1],freshMarAges[x,2]))})
stage.start <- sapply(1:Ncycles,function(x){rep(c(1,2),times=c(freshMarAges[x,1],freshMarAges[x,2]))})
survivorship <- survival <- survivalstage <- matrix(0,nrow=Ncycles,ncol=length(popDyn[1,1,]))
fecundity <- matrix(0,nrow=Ncycles,ncol=length(popDyn[1,1,]))
for(i in 1:Ncycles){
  survival[i,1:length(surv.start[[i]])] <- surv.start[[i]]
  survivalstage[i,1:length(surv.start[[i]])] <- stage.start[[i]]
  survivorship[i,1:length(surv.start[[i]])] <- cumprod(c(1,surv.start[[i]][-length(surv.start[[i]])]))
  fecundity[i,length(surv.start[[i]])] <- 1
}

# leading parameters
SPR0 <- sum(survivorship*fecundity*propAge) #spawner per recruit
CR <- 4 # compensation ratio
N0 <- 100# equilibrium adults
R0 <- N0/SPR0 # equilibrium recruits
alpha.R <- CR/SPR0
beta.R <- -log(alpha.R*SPR0)/(R0*SPR0)
curve(alpha.R*x*exp(beta.R*x),from=0,to=R0*SPR0,xlab="Spawners",ylab="Recruits")

popDyn[1,,] <- R0*survivorship*propAge
Spawners[1] <- sum(sapply(1:Ncycles,function(x){popDyn[1,x,sum(freshMarAges[x,])]}))
Recruits[1] <- alpha.R*Spawners[1]*exp(beta.R*Spawners[1])

RecruitsLC[1,] <- alpha.R*Spawners[1]*exp(beta.R*Spawners[1])*propAge
SpawnersLC[1,] <- sapply(1:Ncycles,function(x){popDyn[1,x,sum(freshMarAges[x,])]})
for(Iyear in 2:Nyears)
{
  # get recruitment from last years' spawners
  Spawners[Iyear] <- max(0,sum(sapply(1:Ncycles,function(x){popDyn[Iyear-1,x,sum(freshMarAges[x,])]})))
  expectedRec <- alpha.R*Spawners[Iyear]*exp(beta.R*Spawners[Iyear])
  Recruits[Iyear] <- rlnorm(1,log(expectedRec),recCV)
  RecruitsLC[Iyear,] <- Recruits[Iyear]*propAge
  SpawnersLC[Iyear,] <- sapply(1:Ncycles,function(x){popDyn[Iyear-1,x,sum(freshMarAges[x,])]})
  
  # calculate time-varying freshwater and marine survival
  frbetaPars <- get_beta(freshSurv,freshSurvCV)
  freshSurvYr[Iyear] <- rbeta(1,frbetaPars$alpha,frbetaPars$beta)
  freshSurvYr[Iyear] <- (freshRho)*freshSurvYr[Iyear-1]+(1-freshRho)*freshSurvYr[Iyear]
  marbetaPars <- get_beta(marSurv,marSurvCV)
  marSurvYr[Iyear] <- rbeta(1,marbetaPars$alpha,marbetaPars$beta)
  marSurvYr[Iyear] <- (marRho)*marSurvYr[Iyear-1]+(1-marRho)*marSurvYr[Iyear]
  for(Ilife in 1:Ncycles)
  {
    popDyn[Iyear,Ilife,1] <- Recruits[Iyear]*propAge[Ilife]
    for(Iage in 2:length(popDyn[Iyear,Ilife,]))
    {
      surv <- ifelse(survivalstage[Ilife,Iage]==1,freshSurvYr[Iyear],
                     ifelse(survivalstage[Ilife,Iage]==2,marSurvYr[Iyear],0))
      popDyn[Iyear,Ilife,Iage] <- popDyn[Iyear-1,Ilife,Iage-1]*surv
    }
  }
}
closureRisk <- sum(ifelse(Spawners>(0.1*N0),0,1))/Nyears

plot(Spawners,Recruits,xlab="Spawners",ylab="Recruits (age 1)")
matplot(SpawnersLC,type="l",xlab="Year",ylab="Spawners",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5)
lines(rowSums(SpawnersLC))
plot(freshSurvYr,type="l",ylab="Freshwater survival (yr)",xlab="Year")
plot(marSurvYr,type="l",ylab="Marine survival (yr)",xlab="Year")
