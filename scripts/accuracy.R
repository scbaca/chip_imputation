# summarize accuracy of imputation for het snps, with array-based genotypes as ground truth
# takes as a commandline argument the *compare.snps.txt files generated scripts/extract_snps.sh
library(stringr)

args=commandArgs(T)
t = read.table(args[1],sep="\t",header=FALSE)

# require that a genotype call was made in the gold-standard vcf
t = t[t[,5]!="./.",]

overlap = function(c, p, chr, start, end) {
        any((c == chr) & (p >= start & p <= end))
}

# filter to include only SNPs overlapping the intervals in the specified bed file
if (length(args)==2) {
	filt = read.table(args[2],sep="\t",header=F)
	filt[,1] = str_replace(filt[,1], "chr", "")

	include=rep(F, nrow(t))
	for (i in 1:nrow(t)) {
		include[i] = overlap(c = t[i, 1], p = t[i,2], chr = filt[,1], start = filt[,2], end = filt[,3])
	}	
	t = subset(t, include)
}

x = t[,5] == "0|1" | t[,5] == "1|0" # find true hets 
y = t[,4] # gene dosage from sequencing data

for ( s in seq(0,1,0.01) ) { # test different cutoffs for gene dosage
keep = abs(y-1) < s
sens = sum(x[keep]) / sum(x)
spec = mean( x[keep] ) 
cat( s , sens , spec , '\n' , sep='\t' )
}
