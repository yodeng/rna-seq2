#******************************************************************part 1 usage ***************************************

library('getopt');
library(ggplot2)
library(reshape2)
require(RColorBrewer)
spec = matrix(c(
  'help'    , 'h', 0, "logical",
  'fpkm'  , 'f', 1, "character",
  'count'  , 'c', 1, "character",
  'o'  , 'o', 1, "character",
  'pre', 'p', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      Rscript  fpkm_saturation_v3.R  --fpkm  genes.fpkm_table  --count genes.count_table  --o result/
      
      Options: 
      --help     -h 	NULL 		get this help
      --fpkm	character 	    the genes.fpkm_table  [forced]
      --count	character 	    the genes.count_table  [forced]
      --o 		character 	    the outfile dirname    [forced]
      \n")
  q(status=1);
}
################################################################################################
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$fpkm) )  { print_usage(spec) }
if ( is.null(opt$count) )  { print_usage(spec) }
if ( is.null(opt$o) )	{ print_usage(spec) }
################################################################

width=140
height=140
size_ratio=width/210

fpkms=read.table(opt$fpkm,sep="\t",header=T,comm="",check=F,row=1)
counts=read.table(opt$count,sep="\t",header=T,comm="",check=F,row=1)
print(head(fpkms))
print(head(counts))
sam=colnames(fpkms)
samples=sam
sample=unlist(lapply(samples, function(x) x[1]))

for (i in 1:length(sample)){

	Total=sum(counts[,i])
	Len=(counts[,i]/(fpkms[,i]*Total))*1e9
	Len[is.na(Len)]=1
	prob=counts[,i]/Total
	count=rep(0,nrow(fpkms))
	names(count)=rownames(fpkms)

	mysample=function(rate){
	  s=sample(rownames(fpkms),size=ceiling(Total*rate),re=T,prob=prob)
	  ct=table(s)
	  count[names(ct)]=ct
	  count/sum(count)/Len*1e9
	}

	d=cbind(sapply(seq(0.1,0.9,0.1),mysample),fpkms[,i])
	#final=d[,10]

	#Number<-apply(d,2,function(x){
		#length(which(abs(x-final)/final<0.15 & final>1 & x>1 ))
		#length(which(x>0.1))
	#})



	 # saturation=data.frame(precent=seq(10,100,10),Number,check.names =F)
	 #p=ggplot(saturation,aes(x=precent,y=Number))+  geom_line(colour="#00a651")+ geom_point(size=1.3,colour="#00a651")+ ylab("Gene Number")+ xlab("Mapped Reads(%)")+ggtitle("Rarefaction Curve")
	# p=p+theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
	# ggsave(filename=paste(opt$o,"/",sample[i],".saturation1.png",sep=""),width=width,height=height,unit="mm",dpi=300)


	#------------------------------------------------------------------------
	d=d[d[,10]>0.1,]


	#--------------------------------------------------------------
	#calculate the absolution difference percent out of final fpkm
	#--------------------------------------------------------------
	ch=apply(d,1,function(x){
	  abs((x-x[10])/x[10])
	})
	ch=t(ch)

	#------------------------------------------------------------------
	#the gene number with  diff percnet <0.15 in ervery interval fpkm
	#-----------------------------------------------------------------
	fpkm=d[,10]
	box<-summary(fpkm)
	print(box)	
	interval=apply(ch,2,function(x){
	  fm0=fpkm[x<0.15]   

	  c(length(which(fm0<=box[2])),
	    length(which(fm0<=box[3]&fm0>box[2])),
	    length(which(fm0<=box[4]&fm0>box[3])),
	    length(which(fm0>box[5]))
	  )
	  
	})

	#-------------------------------------------------------
	#the percent gene out of final value
	#-------------------------------------------------------
	interval=apply(interval,1,function(x){x/x[10]})
	colnames(interval)=c("0-25%","25%-50%","50%-75%","75%-100%")
	data=data.frame(precent=seq(10,100,10),interval,check.names =F)

	data=melt(data,id=1)

	#----------------------------------------------------------------
	#ggplot 
	#---------------------------------------------------------------

	mytheme=theme_get()
	mytheme$axis.title.x$size=14*size_ratio
	mytheme$axis.title.y$size=14*size_ratio
	mytheme$axis.text.x$size=12*size_ratio
	mytheme$axis.text.y$size=12*size_ratio
	mytheme$legend.title$size=11*size_ratio
	mytheme$legend.text$size=9*size_ratio
	theme_set(mytheme)

	p=ggplot(data,aes(x=precent,y=value,colour=variable))+
	  geom_line()+
	  geom_point(size=2*width/210)+
	  ylab("Percent of genes within 15% float of final FPKM")+
	  xlab("Mapped Reads(%)")+
	  scale_colour_manual(values =c("#cc2424","#b276b2","#5da5da","#60bd68"),name="FPKM interval")
	p=p+theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
	ggsave(filename=paste(opt$o,"/",sample[i],".saturation.png",sep=""),width=width,height=height,unit="mm",dpi=300)
}







