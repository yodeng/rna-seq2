library(ggplot2)
library(reshape2)
library(RColorBrewer)

args <- commandArgs(TRUE)

dataFile <- args[1]
outPngNum <- args[2]
samplename <- args[3]

datafile<-read.table(dataFile,header=TRUE,sep="\t")
Sample=datafile$Sample
log2Fpkm=datafile$log2Fpkm
pic<-ggplot(datafile,aes(x=Sample,y=log2Fpkm))+geom_point(position="jitter",col=2,pch=16,cex=1)+stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21,size = 2.5)+ geom_violin(colour='black',fill='gold')+geom_boxplot(width=0.1,col="black",outlier.colour=NA)+theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))
ggsave(pic,file=outPngNum)

