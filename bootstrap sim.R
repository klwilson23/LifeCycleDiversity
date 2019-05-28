source("Some Functions.R")
source("Life cycle function.R")

Ncycles <- 4
Nyears <- 100
Nboots <- 1000

freshMarAges <- matrix(NA,ncol=2,nrow=Ncycles)
freshMarAges[1,] <- c(1,2)
freshMarAges[2,] <- c(1,3)
freshMarAges[3,] <- c(2,3)
freshMarAges[4,] <- c(2,2)

lifeCycleNames <- sapply(1:Ncycles,function(x){paste(freshMarAges[x,],collapse=".")})
dimnames(freshMarAges)=list("Life cycles"=lifeCycleNames,
                            "Stage"=c("Freshwater","Marine"))

propAge1 <- rep(1/Ncycles,Ncycles)
propAge2 <- c(0.1,0.1,0.1,0.7)
propAge3 <- c(0.1,0.1,0.7,0.1)
propAge4 <- c(0.1,0.7,0.1,0.1)
propAge5 <- c(0.7,0.1,0.1,0.1)
propAge <- rbind(propAge1,propAge2,propAge3,propAge4,propAge5)

recCV <- 0.35

freshSurvMn <- 0.075
freshSurvCV <- 0.25
freshRho <- 0.8

marSurvMn <- 0.075
marSurvCV <- 0.25
marRho <- 0.8

CR <- 6 # compensation ratio
N0 <- 1e6/0.2 # equilibrium adults
propRisk <- 0.2 # fishery closes if population falls below propRisk*N0 

risk <- c()
for(i in 1:nrow(propAge))
{
  boot <- replicate(Nboots,LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[i,],freshMarAges=freshMarAges,recCV=recCV,Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk))
  risk[i] <- mean(unlist(boot[which(dimnames(boot)[[1]]=="closureRisk"),]))
}

risk/risk[1]

example <- LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[1,],freshMarAges=freshMarAges,recCV=recCV,Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk)

plot(example$Spawners,example$Recruits,xlab="Spawners",ylab="Recruits")
matplot(example$SpawnersLC,type="l",xlab="Year",ylab="Spawners",col=c("orange","dodgerblue","grey50","blue"),lty=3,lwd=1.5)
lines(rowSums(example$SpawnersLC))
plot(example$freshSurvYr,type="l",ylab="Freshwater survival",xlab="Year")
plot(example$marSurvYr,type="l",ylab="Marine survival",xlab="Year")
