txtSize <- 0.5
Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}


jpeg("comparing time periods.jpeg",width=7,height=4,units="in",res=800)
escape1 <- 0.88e6
escape2 <- 0.88e6
cv1 <- 1.005
cv1_log <- sqrt(log(1+(cv1*escape1)^2/(escape1^2)))
mu1_log <- log((escape1^2)/sqrt(escape1^2+(cv1*escape1)^2))
cv2 <- 0.677
cv2_log <- sqrt(log(1+(((cv2*escape2)^2)/(escape2^2))))
thresh <- 1.05e6
escape1.boot <- pmax(0,dlnorm(1e6,meanlog=log(escape1)-cv1_log^2/2,sdlog=cv1_log))
escape2.boot <- pmax(0,rlnorm(1e6,meanlog=log(escape2)-cv2_log^2/2,sdlog=cv2_log))
escape3 <- 1.3e6
escape4 <- 1.3e6
cv3 <- 1.005
cv3_log <- log(1+(cv3*escape3)^2/(escape3^2))
cv4 <- 0.677
cv4_log <- log(1+(cv4*escape4)^2/(escape4^2))
thresh <- 1.05e6
escape3.boot <- pmax(0,rlnorm(1e6,meanlog=log(escape3)-cv3_log^2/2,sdlog=cv3_log))
escape4.boot <- pmax(0,rlnorm(1e6,meanlog=log(escape4)-cv4_log^2/2,sdlog=cv4_log))


brk.seq <- seq(min(escape1.boot,escape2.boot,escape3.boot,escape4.boot)-1,max(escape1.boot,escape2.boot,escape3.boot,escape4.boot)+1,length=100)

layout(matrix(1:2,ncol=2))
par(mar=c(3,2,1,0),mai=c(0.6,0.6,0.1,0.05))
hist(escape1.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,26e-7),main=paste("Mean abundance = ",escape1,sep=""),xlab="",cex.axis=txtSize,cex.lab=txtSize,cex.main=txtSize,xlim=c(0,max(brk.seq)),ylab="")
mtext("Abundance",side=1,line=2,cex=txtSize)
mtext("Density",side=2,line=2,cex=txtSize)
hist(escape2.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos1 <- round(100*sum(escape1.boot<=thresh)/length(escape1.boot),0)
perc_clos2 <- round(100*sum(escape2.boot<=thresh)/length(escape2.boot),0)
perc_extinc1 <- round(100*sum(escape1.boot<=0)/length(escape1.boot),0)
perc_extinc2 <- round(100*sum(escape2.boot<=0)/length(escape2.boot),0)
Corner_text(paste(perc_clos1,"% high variance closures","\n",perc_clos2,"% low variance closures","\n\n",perc_extinc1,"% high variance extirpation","\n",perc_extinc2,"% low variance extirpation",sep=""),"topright",cex=txtSize)
abline(v=thresh,lwd=3*txtSize,lty=3,col="red")

par(mar=c(0,2,1,3),mai=c(0.6,0.05,0.1,0.6))
hist(escape3.boot,breaks=brk.seq,col=adjustcolor("dodgerblue",0.5),freq=FALSE,ylim=c(0,26e-7),main=paste("Mean abundance = ",escape3,sep=""),xlab="",cex.axis=txtSize,cex.lab=txtSize,cex.main=txtSize,xlim=c(0,max(brk.seq)),yaxt="n")
mtext("Abundance",side=1,line=2,cex=txtSize)

hist(escape4.boot,breaks=brk.seq,col=adjustcolor("blue",0.5),add=TRUE,freq=FALSE)
perc_clos3 <- round(100*sum(escape3.boot<=thresh)/length(escape3.boot),0)
perc_clos4 <- round(100*sum(escape4.boot<=thresh)/length(escape4.boot),0)
perc_extinc3 <- round(100*sum(escape3.boot<=0)/length(escape3.boot),0)
perc_extinc4 <- round(100*sum(escape4.boot<=0)/length(escape4.boot),0)
Corner_text(paste(perc_clos3,"% high variance closures","\n",perc_clos4,"% low variance closures","\n\n",perc_extinc3,"% high variance extirpation","\n",perc_extinc4,"% low variance extirpation",sep=""),"topright",cex=txtSize)
legend("right",c("High variance","Low variance"),pch=22,pt.bg=c("dodgerblue","blue"),bty="n",title="Scenarios",cex=txtSize)
abline(v=thresh,lwd=3*txtSize,lty=3,col="red")
dev.off()
