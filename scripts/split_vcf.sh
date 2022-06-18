# split a vcf with many samples into many vcfs with one sample each
# this is used to split the "reference" calls from Blueprint blood samples
# eg: sh scripts/split_vcf.sh ../../../EGA/blood/EGAF00001331960/All_chr.BPWP10_13_12_15.vcf.gz blueprint_gt

file=$1
outdir=$2
mkdir -p $outdir
for sample in `bcftools query -l $file`; do
	outfile=${outdir}/${sample}.vcf.gz
	if [ -f "$outfile" ]; then
		echo skipping $outfile, which already exists
	else
		echo creating file $outfile
		bcftools view -Oz -s $sample -o $outfile $file
	fi
done
