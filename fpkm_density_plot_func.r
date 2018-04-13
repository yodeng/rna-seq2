
# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_density_plot_func.r  genes.fpkm_table  outdir")
	print("1) genes.fpkm_table: the fpkm file for RNA-SEQ")
	print("2) outdir: the dir for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
	print(args)
	usage()
	stop("the length of args != 2")
}


# load library
require(edgeR)
require(ggplot2)


fpkms=read.table(args[1],sep="\t",header=T,comm="",check=F,row=1)
sam=colnames(fpkms)
samples=sam
sam_name=unlist(lapply(samples, function(x) x[1]))

print(head(fpkms))
all<-NULL
all_sam<-NULL
# iter plot fpkm density
for( i in 1:length(sam_name) ){
	keep <- fpkms[,i] > 0
	r <-fpkms[keep,i]
	print(head(r))
	log10fpkm <- data.frame(log10fpkm=log10(r))

	# update all
	all <- c(all, log10(r))
	all_sam <- c(all_sam, rep(sam_name[i], length(r)))

	# plot
	m <- ggplot(log10fpkm, aes(x=log10fpkm))
	p <- m + geom_density(fill="#00a651", size=1, colour="#00a651") + xlab("log10(FPKM)")
	p <- p + theme(axis.title.x = element_text(face="bold", size=14),
		axis.text.x  = element_text(face="bold", size=12),
		axis.title.y = element_text(face="bold", size=14),
		axis.text.y  = element_text(face="bold", size=12) )
	p <- p + theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=12) )
	p <- p+theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
	# output plot into file
	file <- paste(args[2], "/", sam_name[i], ".fpkm_density.png", sep="")
	png(filename=file, height = 3000, width = 3000, res = 500, units = "px")
	print(p)
	dev.off()
}


# plot fpkm density for all
# create data.frame
log10fpkm <- data.frame(log10fpkm=all, sample=all_sam)
Sample <- factor(all_sam)
# plot
m <- ggplot(log10fpkm, aes(x=log10fpkm))
p <- m + geom_density(aes(fill=Sample, colour=Sample),alpha = 0.3) + xlab("log10(fpkm)")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
p<-p+theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())

if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
# output
file <- paste(args[2], "/all", ".fpkm_density.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()


# plot fpkm box for all
# plot
m <- ggplot(log10fpkm, aes(factor(sample), log10fpkm))
p <- m + geom_boxplot(aes(fill=Sample)) + xlab("Sample")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
 p<-p+theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
# output
file <- paste(args[2], "/all", ".fpkm_box.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()







