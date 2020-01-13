Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

escape1 <- 0.88e6
escape2 <- 0.88e6
cv1 <- 1.005
cv2 <- 0.677
thresh <- 1.05e6
escape1.boot <- pmax(0,rnorm(1e6,mean=escape1,sd=cv1*escape1))
escape2.boot <- pmax(0,rnorm(1e6,mean=escape2,sd=cv2*escape2))

brk.seq <- seq(min(escape1.boot,escape2.boot)-1,max(escape2.boot,escape1.boot)+1,length=100)

layout(matrix(1:2,ncol=2))
hist(escape1.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,33e-7),main=paste("Escapement for high- versus low-variance","\n mean N = ",escape1,sep=""),xlab="Escapement")
hist(escape2.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos1 <- round(100*sum(escape1.boot<=thresh)/length(escape1.boot),0)
perc_clos2 <- round(100*sum(escape2.boot<=thresh)/length(escape2.boot),0)
Corner_text(paste(perc_clos1,"% high variance closures","\n",perc_clos2,"% low variance closures",sep=""),"topright")
legend("right",c("Higher variance, same abundance","Lower variance, same abundance"),pch=22,pt.bg=c("dodgerblue","blue"),bty="n",title="Scenarios")
abline(v=thresh,lwd=3,lty=3,col="red")

escape1 <- 1.3e6
escape2 <- 1.3e6
cv1 <- 1.005
cv2 <- 0.677
thresh <- 1.05e6
escape1.boot <- pmax(0,rnorm(1e6,mean=escape1,sd=cv1*escape1))
escape2.boot <- pmax(0,rnorm(1e6,mean=escape2,sd=cv2*escape2))

brk.seq <- seq(min(escape1.boot,escape2.boot)-1,max(escape2.boot,escape1.boot)+1,length=100)

hist(escape1.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,33e-7),main=paste("Escapement for high- versus low-variance","\n mean N = ",escape1,sep=""),xlab="Escapement")
hist(escape2.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos1 <- round(100*sum(escape1.boot<=thresh)/length(escape1.boot),0)
perc_clos2 <- round(100*sum(escape2.boot<=thresh)/length(escape2.boot),0)
Corner_text(paste(perc_clos1,"% high variance closures","\n",perc_clos2,"% low variance closures",sep=""),"topright")
legend("right",c("Higher variance, same abundance","Lower variance, same abundance"),pch=22,pt.bg=c("dodgerblue","blue"),bty="n",title="Scenarios")
abline(v=thresh,lwd=3,lty=3,col="red")