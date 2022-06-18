#plot ancestry principal compoments: 
#run as, eg: Rscript plot.R output/CO.PC1.profile output/CO.PC2.profile CO.PC1.PC2
#This is not integrated into the workflow currently

library(stringr)

args=commandArgs(T)
t1 = read.table(args[1],header=TRUE)
t2 = read.table(args[2],header=TRUE)
plot.name = args[3]

mat = data.frame("pc1"=t1$SCORE, "pc2"=t2$SCORE)
mat$sample=t1$IID

pc1.mean=mean(mat$pc1)
pc1.sd=sd(mat$pc1)
pc2.mean=mean(mat$pc2)
pc2.sd=sd(mat$pc2)

mat$outlier = mat$pc1 > pc1.mean + 2 * pc1.sd | mat$pc1 < pc1.mean - 2 * pc1.sd |
  mat$pc2 > pc2.mean + 2 * pc2.sd | mat$pc2 < pc2.mean - 2 * pc2.sd

#add self-reported ancestry annotations:
#ann=read.table("annotations.txt", head=T, as.is=T)
#mat=merge(mat, ann, by="sample", all.x=T)
#mat$Ancestry[is.na(mat$Ancestry)]="European or unspecified"
#mat$Ancestry=as.factor(mat$Ancestry)

library(ggplot2)
ggplot(mat,aes(x=pc1,y=pc2)) + geom_point() + theme_classic() + geom_text(data=mat[mat$outlier,], aes(label=sample), show_guide=F) + xlab("PC1") + ylab("PC2") 

ggsave(plot.name, height=4, width=6)


#print outliers to file:
outliers = mat$sample[mat$outlier]

write.table(outliers, str_replace(plot.name, '.pdf', '.outliers'), row.names=F, col.names=F, quote=F)
