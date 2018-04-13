library(ggplot2)
library(reshape2)
library(RColorBrewer)

args <- commandArgs(TRUE)

dataFile <- args[1]
outPngNum <- args[2]
outPngNumlog <- args[3]
outPngLength <- args[4]
outPngLengthlog <- args[5]

datafile<-read.table(dataFile,header=TRUE,sep="\t")
Type=datafile$Type
log2Transcript_exon_num=datafile$log2Transcript_exon_num
log2Transcript_length=datafile$log2Transcript_length
Transcript_exon_num=datafile$Transcript_exon_num
Transcript_length=datafile$Transcript_length

pic<-ggplot(datafile,aes(x=Type,y=Transcript_exon_num))+geom_point(position="jitter",col=2,pch=16,cex=1)+geom_boxplot(col="black",pch=16,cex=1)+theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))
ggsave(pic,file=outPngNum)
pic<-ggplot(datafile,aes(x=Type,y=log2Transcript_exon_num))+geom_point(position="jitter",col=2,pch=16,cex=1)+geom_boxplot(col="black",pch=16,cex=1)+theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))
ggsave(pic,file=outPngNumlog)


pic1<-ggplot(datafile,aes(x=Type,y=Transcript_length))+geom_point(position="jitter",col=2,pch=16,cex=1)+geom_boxplot(col="black",pch=16,cex=1)+theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))
ggsave(pic1,file=outPngLength)
pic1<-ggplot(datafile,aes(x=Type,y=log2Transcript_length))+geom_point(position="jitter",col=2,pch=16,cex=1)+geom_boxplot(col="black",pch=16,cex=1)+theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))
ggsave(pic1,file=outPngLengthlog)
