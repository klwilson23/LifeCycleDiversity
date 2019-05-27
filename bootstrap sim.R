source("Some Functions.R")
source("Life cycle function.R")

Ncycles <- 4
Nyears <- 100
Nboots <- 1000

freshMarAges <- matrix(NA,ncol=2,nrow=Ncycles)
freshMarAges[1,] <- c(2,3)
freshMarAges[2,] <- c(3,2)
freshMarAges[3,] <- c(2,4)
freshMarAges[4,] <- c(3,3)

lifeCycleNames <- sapply(1:Ncycles,function(x){paste(freshMarAges[x,],collapse=".")})
dimnames(freshMarAges)=list("Life cycles"=lifeCycleNames,
                            "Stage"=c("Freshwater","Marine"))

propAge1 <- rep(1/Ncycles,Ncycles)
propAge2 <- c(0.1,0.1,0.1,0.7)
propAge <- rbind(propAge1,propAge2)

recCV <- 0.35

freshSurvMn <- 0.1
freshSurv <- exp(--log(freshSurvMn)/sum((propAge*freshMarAges[,1])))
freshSurvCV <- 0.25
freshRho <- 0.6

marSurvMn <- 0.075
marSurv <- exp(--log(marSurvMn)/sum((propAge*freshMarAges[,2])))
marSurvCV <- 0.25
marRho <- 0.6

CR <- 4 # compensation ratio
N0 <- 100 # equilibrium adults
propRisk <- 0.25 # fishery closes if population falls below propRisk*N0 

risk <- c()
for(i in 1:nrow(propAge))
{
  boot <- replicate(Nboots,LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[i,],freshMarAges=freshMarAges,recCV=recCV,Ncycles=Ncycles,freshSurv=freshSurv,marSurv=marSurv,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk))
  risk[i] <- mean(unlist(boot[which(dimnames(boot)[[1]]=="closureRisk"),]))
}

risk[2]/risk[1]
