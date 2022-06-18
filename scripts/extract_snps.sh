# merge SNP calls from imputation and arrays for a given sample

imp=$1 # vcf file with imputed genotypes
array=$2 # vcf file with array-based genotypes
out=$3 # output folder
arrayid=$4 # individual id in the snp array vcf
name=$5 # individaul id

mkdir -p $out

echo extracting imputed genotypes for $name
bcftools view $imp -Oz -s $name > ${out}/${name}.imp.tmp.vcf.gz

echo extracting array-based genotypes for $name with SNP array ID $arrayid
bcftools view $array -Oz -s $arrayid > ${out}/${name}.array.tmp.vcf.gz

echo indexing vcf files
bcftools index ${out}/${name}.imp.tmp.vcf.gz
bcftools index ${out}/${name}.array.tmp.vcf.gz

echo merging imputed and reference genotype files
bcftools merge -m none --force-samples ${out}/${name}.imp.tmp.vcf.gz ${out}/${name}.array.tmp.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT\t%DS\t]\n' |  awk 'BEGIN{OFS="\t"} $4!="." {print $0}' > ${out}/${name}.compare.snps.txt

rm ${out}/${name}.imp.tmp.vcf.gz ${out}/${name}.array.tmp.vcf.gz

