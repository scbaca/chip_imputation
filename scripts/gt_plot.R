library(ggplot2)
library(stringr)
library(reshape2)

args = commandArgs(trailingOnly = TRUE)
gt.file = args[1]
out.dir = args[2]
metasheet = args[3]

meta = read.table(metasheet, header=T, sep=",")
mat = read.table(gt.file,header = F, colClasses = c("numeric","numeric","character","character"))

colnames(mat) = c("error", "n.snps", "sample1", "sample2")

mat$sample1 = str_replace(mat$sample1,".*/","")
mat$sample1 = str_replace(mat$sample1,"_unique.sorted.dedup.bam","")

mat$sample2 = str_replace(mat$sample2,".*/","")
mat$sample2 = str_replace(mat$sample2,"_unique.sorted.dedup.bam","")

#prepend individual names to each sample name
meta$bam = str_replace(meta$bam,".*/","")
meta$bam = str_replace(meta$bam,"_unique.sorted.dedup.bam","")

colnames(meta)[1] = "sample1"
mat = merge(mat, meta, by="sample1")
mat$sample1 = paste(mat$individual, mat$sample1, sep = "_")
mat = mat[,1:ncol(mat)-1] #remove the "individual" row

colnames(meta)[1] = "sample2"
mat = merge(mat, meta, by="sample2")
mat$sample2 = paste(mat$individual, mat$sample2, sep = "_")
mat = mat[,1:ncol(mat)-1] #remove the "individual" row

#tmp - tidy up names for blueprint samples
mat$sample1 = mat$sample1 %>% str_replace(".bwa.*","") %>% str_replace("\\.", "_") %>% str_replace ("NS.*_", "")
mat$sample2 = mat$sample2 %>% str_replace(".bwa.*","") %>% str_replace("\\.", "_") %>% str_replace ("NS.*_", "")

#add a self to self comparisons:
all.samples = unique(c(mat$sample1,mat$sample2))


n = length(all.samples) 
append = data.frame("error" = rep(0,n), "n.snps" = rep(0,n), "sample1" = all.samples, "sample2" = all.samples)
mat = rbind(mat,append)

#include all samples as sample 1 and as sample 2 (ie make a square rather than triangle)
#this fixes issues with some of the comparisons swapping samples 1 and 2 inappropriately in 
#gtcheck output

mat = unique(rbind(mat, data.frame("error" = mat$error, "n.snps" = mat$n.snps, 
	"sample1" = mat$sample2, "sample2" = mat$sample1)))

#convert to matrix for clustering
mat.2d = acast(mat,sample1 ~ sample2, value.var = "error")


#cluster
c = hclust(as.dist(mat.2d))
pdf(file.path(out.dir, "gt.clust.pdf"), height = n * 0.25+3, width = 12)
par(pin=c(3, n*0.25)) #tmp
plot(as.dendrogram(c), cex = 1, main = "sample clustering based on comparison of genotypes at homozygous SNPs", horiz = T, xlab = "distance")
dev.off()

mat$sample1 = factor(mat$sample1, levels = c$labels[c$order])
mat$sample2 = factor(mat$sample2, levels = c$labels[c$order])

# sample-sample correlation heatmap
max.err = max(mat$error)
min.err = min(mat$error[mat$error > 0])

ggplot(mat) + geom_tile(aes(x = sample1, y = sample2, fill = error), color = "white") + 
  scale_fill_gradient(limits = c(min.err, max.err), low = "blue", high = "gold") + 
  theme_minimal(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(out.dir, "gt.heatmap.pdf"), width=n * 0.25, height = n * 0.25, limitsize = FALSE)
