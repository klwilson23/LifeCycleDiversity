source("Some Functions.R")
source("Life cycle function.R")

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

scenarios <- array(NA,dim=c(nrow(propAge),
                           length(recCV),
                           length(probGood),
                           length(probBad)),
                  dimnames=list("age structure"=row.names(propAge),
                                "recCV"=names(recCV),
                                "probGood"=probGood,
                                "probBad"=probBad))
nScenarios <- length(scenarios) # how many scenarios are we simulating

counter <- 0
results <- as.data.frame(expand.grid(dimnames(scenarios)))
results$probGood <- as.numeric(as.character(results$probGood))
results$probBad <- as.numeric(as.character(results$probBad))
results <- data.frame(results,"mn_risk"=rep(NA,nrow(results)),
                      "LI_risk"=rep(NA,nrow(results)),
                      "UI_risk"=rep(NA,nrow(results)))

for(i in 1:nrow(results))
{
  boot <- replicate(Nboots,LifeCycle(Nyears=Nyears,CR=CR,N0=N0,propAge=propAge[results$age.structure[i],],freshMarAges=freshMarAges,recCV=recCV[results$recCV[i]],Ncycles=Ncycles,freshSurvMn=freshSurvMn,marSurvMn=marSurvMn,freshSurvCV=freshSurvCV,marSurvCV=marSurvCV,freshRho=freshRho,marRho=marRho,lifeCycleNames=lifeCycleNames,propRisk=propRisk,marSurv_Scen="g-g",probGood=results$probGood[i],probBad=results$probBad[i],goodSurv=goodSurv,badSurv=badSurv,startProb=startProb))
  results$mn_risk[i] <- mean(unlist(boot[which(dimnames(boot)[[1]]=="closureRisk"),]))
  results$LI_risk[i] <- quantile(unlist(boot[which(dimnames(boot)[[1]]=="closureRisk"),]),probs=0.025)
  results$UI_risk[i] <- quantile(unlist(boot[which(dimnames(boot)[[1]]=="closureRisk"),]),probs=0.975)
}

results$portfolio <- 1/(results$mn_risk/results[results$age.structure=="Recent" & results$recCV=="Historical" & results$probGood==0.8 & results$probBad==0.8,"mn_risk"])

saveRDS(results,file="results_boot.rds")
