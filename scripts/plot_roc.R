#plot ROC curve for detection of het SNPs from imputed data, with array-based genotypes as ground truth

library(stringr)
library(scales)

args = commandArgs(trailingOnly = TRUE)
list = args[1:(length(args)-1)]
plotfile = args[length(args)]
statfile = str_replace(plotfile, ".pdf", ".stats")

names=str_replace(list,".*/","")
names=str_replace(names,".roc.txt","")

pdf(plotfile,width=5,height=5)

stats = data.frame(name=names, sens80=rep(0,length(list)), sens90=rep(0,length(list))) #sensitivity at specificity of 0.8

plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),main="ROC for het SNP detection",xlab="1-Specificity",ylab="Sensitivity")
for ( i in 1:length(list) ) {

	line.col=alpha(rgb(0,0,0), 0.1)

	acc = read.table(list[i], as.is=T )
	lines( 1 - acc[,3] , acc[,2] , col=line.col)
	
	stats$sens80[i] = acc[max(which(acc[,3] > 0.8)),2]
	stats$sens90[i] = acc[max(which(acc[,3] > 0.9)),2]
}

# replace any NAs with 0s, for a conservative estimate of sens80 and sens90
stats$sens80[is.na(stats$sens80)] = 0
stats$sens90[is.na(stats$sens90)] = 0

abline(h=mean(stats$sens90),col="red",lty=2)
box()
dev.off()

write.table(stats, file=statfile, row.names=F, col.names=T, quote=F)
