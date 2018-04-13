args <- commandArgs(TRUE)
dataFile<-args[1]
outputPngFile<-args[2]
logarg<-""
if(is.element(args[3], c("x", "y", "xy", "yx"))) {
	logarg<-args[3]
}

png(outputPngFile,width =1480, height = 1480)
par(mar=c(7,6,6,4))
datafile<-read.table(dataFile,header=FALSE,sep="\t")
header<-subset(datafile[1,])
xlab<-header$V2
ylab<-header$V3
main<-header$V4
datafile<-datafile[-1,]
datafile<-datafile[order(factor(datafile$V4,levels = c("","up","down"))),]
datafile<-as.matrix(datafile)
col1=datafile[,4]
col1[col1=="up"]="#cc2424"
col1[col1=="down"]="#00a651"
col1[col1==""]="#000000"
plot(datafile[,2],datafile[,3],pch=20,cex=3,xlab=xlab,
     ylab=ylab,main=main,col=as.vector(col1),cex.axis=2.25,cex.lab=2.7,cex.main=4,log=logarg)
legend("topleft",pch=20,col=c("#cc2424","#00a651","#000000"),legend=c("Up regulated genes","Down regulated genes","Non-regulated genes"),
       bty="n",xpd=TRUE,cex=2)
dev.off()
