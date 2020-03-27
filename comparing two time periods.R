txtSize <- 0.5
Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}


jpeg("comparing time periods.jpeg",width=7,height=4,units="in",res=800)
escape1 <- 0.88e6
escape2 <- 0.88e6
cv1 <- 1.005
cv2 <- 0.677
thresh <- 1.05e6
escape1.boot <- pmax(0,rnorm(1e6,mean=escape1,sd=cv1*escape1))
escape2.boot <- pmax(0,rnorm(1e6,mean=escape2,sd=cv2*escape2))

brk.seq <- seq(min(escape1.boot,escape2.boot)-1,max(escape2.boot,escape1.boot)+1,length=100)

layout(matrix(1:2,ncol=2))
par(mar=c(3,2,1,0),mai=c(0.6,0.6,0.1,0.05))
hist(escape1.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,35e-7),main=paste("Mean abundance = ",escape1,sep=""),xlab="",cex.axis=txtSize,cex.lab=txtSize,cex.main=txtSize,xlim=c(0,6e6),ylab="")
mtext("Abundance",side=1,line=2,cex=txtSize)
mtext("Density",side=2,line=2,cex=txtSize)
hist(escape2.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos1 <- round(100*sum(escape1.boot<=thresh)/length(escape1.boot),0)
perc_clos2 <- round(100*sum(escape2.boot<=thresh)/length(escape2.boot),0)
Corner_text(paste(perc_clos1,"% high variance closures","\n",perc_clos2,"% low variance closures",sep=""),"topright",cex=txtSize)
abline(v=thresh,lwd=3*txtSize,lty=3,col="red")

escape1 <- 1.3e6
escape2 <- 1.3e6
cv1 <- 1.005
cv2 <- 0.677
thresh <- 1.05e6
escape1.boot <- pmax(0,rnorm(1e6,mean=escape1,sd=cv1*escape1))
escape2.boot <- pmax(0,rnorm(1e6,mean=escape2,sd=cv2*escape2))

brk.seq <- seq(min(escape1.boot,escape2.boot)-1,max(escape2.boot,escape1.boot)+1,length=100)
par(mar=c(0,2,1,3),mai=c(0.6,0.05,0.1,0.6))
hist(escape1.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,35e-7),main=paste("Mean abundance = ",escape1,sep=""),xlab="",cex.axis=txtSize,cex.lab=txtSize,cex.main=txtSize,xlim=c(0,6e6),yaxt="n")
mtext("Abundance",side=1,line=2,cex=txtSize)

hist(escape2.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos1 <- round(100*sum(escape1.boot<=thresh)/length(escape1.boot),0)
perc_clos2 <- round(100*sum(escape2.boot<=thresh)/length(escape2.boot),0)
Corner_text(paste(perc_clos1,"% high variance closures","\n",perc_clos2,"% low variance closures",sep=""),"topright",cex=txtSize)
legend("right",c("Higher variance, same abundance","Lower variance, same abundance"),pch=22,pt.bg=c("dodgerblue","blue"),bty="n",title="Scenarios",cex=txtSize)
abline(v=thresh,lwd=3*txtSize,lty=3,col="red")
dev.off()