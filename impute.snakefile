configfile: "config.yaml"

import pandas as pd
import yaml
import os

meta = pd.read_table(config['metasheet'], sep=',', comment='#')
if config['SNP_array_list'] != '':
	array_list = pd.read_table(config['SNP_array_list'], sep=',', comment='#')

def get_targets(wildcards):
	ls = []
	for ind in meta.individual.unique():
		ls.append('analysis/merged/%s/%s.merged.nochr.bam' % (ind, ind))
		ls.append('analysis/merged/%s/%s.merged.nochr.bam.bai' % (ind, ind))
	ls.append('analysis/stitch/bamlist.txt')
	ls.append('analysis/stitch/indlist.txt')
	ls.append('analysis/genotype/plots/gt.heatmap.pdf')
	ls.append('analysis/stitch/merged/all.genotypes.vcf.gz')
	ls.append('analysis/stitch/merged/all.genotypes.AF.0.05.vcf.gz')
	if config['SNP_array_list'] != '':
		ls.append('analysis/roc/plots/roc.pdf')
	return ls		


def get_bams(wildcards):
	return meta.bam[meta.individual==wildcards.ind]

def get_merged_bams(wildcards):
	ls = []
	for ind in meta.individual.unique():
		ls.append('analysis/merged/%s/%s.merged.nochr.bam' % (ind, ind))
	return ls

def get_merged_vcfs(wildcards):
	ls = []
	for ind in meta.individual.unique():
		ls.append('analysis/genotype/%s/%s.merged.vcf.gz' % (ind, ind))
	return ls

def get_snp_array_vcf(wildcards):
	return array_list.snp_array_vcf[array_list.individual==wildcards.ind]

def get_snp_array_id(wildcards):
	return array_list.snp_array_id[array_list.individual==wildcards.ind].item()

def get_roc_targets(wildcards):
        ls = []
        for ind in array_list.individual:
                ls.append('analysis/roc/%s/%s.roc.txt' % (ind, ind))
        return ls

# remove "chr" if present, leaving just chromosome number
rule remove_chr_and_index:
	input:
		'analysis/merged/{ind}/{ind}.merged.bam'
	output:
		bam = 'analysis/merged/{ind}/{ind}.merged.nochr.bam',
		index = 'analysis/merged/{ind}/{ind}.merged.nochr.bam.bai'
	shell:
		''' samtools view -h {input} | sed "s/chr//g" | samtools view -bS > {output.bam}
		samtools index {output.bam} '''

# merge all bam files for a given individual
rule merge_bams:
	input:
		get_bams
	output:
		temp('analysis/merged/{ind}/{ind}.merged.bam')
	shell:
#		'samtools merge {output} {input}'
		' samtools merge - {input} | samtools addreplacerg -r "@RG\\tID:{wildcards.ind}" - > {output} '

# get genotypes from bam files for qc purposes
rule genotype_from_bam:
	input:
		get_bams
	output:
		'analysis/genotype/{ind}/{ind}.merged.vcf.gz'
	params:
		snps=config['qc_snps'],
		genome=config['genome']
	shell:
		''' mkdir -p analysis/genotype/{wildcards.ind}
		printf '' > analysis/genotype/{wildcards.ind}/merge.list
		printf '' > analysis/genotype/{wildcards.ind}/bams.list
		for f in {input}; do
			name=`echo $f | sed 's/.*\///' | sed 's/.bam//'`
			echo getting genotypes from $name for individual {wildcards.ind}
			samtools index $f
			set +o pipefail; bcftools mpileup -Ou -R {params.snps} -f {params.genome} $f | bcftools call -c | head -n 100000 | bgzip > analysis/genotype/{wildcards.ind}/$name.vcf.gz
			bcftools index analysis/genotype/{wildcards.ind}/$name.vcf.gz
			ls analysis/genotype/{wildcards.ind}/$name.vcf.gz >> analysis/genotype/{wildcards.ind}/merge.list
			ls $f >> analysis/genotype/{wildcards.ind}/bams.list
		done 
		bcftools merge -Oz --force-samples -l analysis/genotype/{wildcards.ind}/merge.list > {output} 
		echo {wildcards.ind} > {output}.rename.tmp
#		bcftools reheader -s {output}.rename.tmp -o {output}.tmp.renamed {output}
		bcftools reheader -s analysis/genotype/{wildcards.ind}/bams.list -o {output}.tmp.renamed {output}
		mv {output}.tmp.renamed {output}
		rm {output}.rename.tmp
		bcftools index {output} '''

# cluster samples by genotype from bam files to make sure there are no sample/individual mixups
rule gt_cluster_plot:
	input:
		vcfs=get_merged_vcfs,
		meta=config['metasheet']
	output:
		'analysis/genotype/plots/gt.heatmap.pdf'
	shell:
		''' bcftools merge -Oz {input.vcfs} > analysis/genotype/plots/merge.all.vcf.gz 
		bcftools gtcheck -H analysis/genotype/plots/merge.all.vcf.gz  | grep $'^ER' | cut -f 2-5 > analysis/genotype/plots/gt.mat
		Rscript scripts/gt_plot.R analysis/genotype/plots/gt.mat analysis/genotype/plots {input.meta} '''

# create list of bam files and names for use in stitch workflow
rule list_bams:
	input:
		bams = get_merged_bams,
# TODO: fix this:		indexes = [b + '.bai' for b in get_merged_bams] #included so indexes are generated before run_stitch rule
	output:
		bamlist = 'analysis/stitch/bamlist.txt',
		indlist = 'analysis/stitch/indlist.txt'
	shell:
		''' ls {input.bams} > {output.bamlist}
		sed 's/.*\///' {output.bamlist} | sed 's/.merged.nochr.bam//' > {output.indlist} ''' 

# impute genotypes with stitch
rule run_stitch:
	input:
		bamlist = 'analysis/stitch/bamlist.txt',
		indlist = 'analysis/stitch/indlist.txt'
	output:
		vcf = 'analysis/stitch/output/stitch.{chr}.{start}.{end}.vcf.gz',
		index = 'analysis/stitch/output/stitch.{chr}.{start}.{end}.vcf.gz.csi',
	params:
		ref_fasta = config['genome'],
		k = config['k'],
		nGen = config['nGen'],
		niterations = config['niterations'],
		buffer = config['buffer'],
		imputation_path = config['imputation_path'],
	shell:
		''' R -e 'library("STITCH"); STITCH(outputdir = "analysis/stitch/output", \
                  bamlist="{input.bamlist}", \
                  reference="{params.ref_fasta}", \
                  posfile="{params.imputation_path}/posfiles/HRC.r1-1.GRCh37.wgs.mac5.sites.{wildcards.chr}.EUR_AF01.posfile", \
                  K={params.k}, \
                  nGen={params.nGen}, \
                  nCores=1, \
                  niterations={params.niterations}, \
                  regionStart={wildcards.start}, \
                  regionEnd={wildcards.end}, \
                  chr={wildcards.chr}, \
                  buffer={params.buffer}, \
                  reference_legend_file="{params.imputation_path}/legend/1000GP_Phase3_chr{wildcards.chr}.legend.gz", \
                  reference_haplotype_file="{params.imputation_path}/hapmap/1000GP_Phase3_chr{wildcards.chr}.hap.gz", \
                  reference_sample_file="{params.imputation_path}/1000GP_Phase3.sample", \
                  reference_populations=c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","POP","PUR","STU","TSI","YRI"), \
                  sampleNames_file="{input.indlist}", \
                  method="diploid", \
                  regenerateInput=TRUE, \
                  inputBundleBlockSize=100, \
                  shuffleHaplotypeIterations=NA, \
                  refillIterations=NA, \
                  output_format="bgvcf", \
                  B_bit_prob=8)'
		bcftools index {output.vcf}
                '''

# merge vcf files generated by stitch
rule merge_stitch:
	input:
		vcfs = config['stitch_targets'],
		indexes = [s + '.csi' for s in config['stitch_targets']]
	output:
		'analysis/stitch/merged/all.genotypes.vcf.gz'
	shell:
		''' ls {input.vcfs} | sort -k2,2n -k3,3n -t '.' > analysis/stitch/merged/merge.list
		bcftools concat -n -f analysis/stitch/merged/merge.list -Oz -o {output} 
		bcftools index {output}'''

# create a vcf file limited to MAF > 0.05 in HRC data
rule subset_vcf_by_MAF:
	input:
		'analysis/stitch/merged/all.genotypes.vcf.gz'
	output:
		'analysis/stitch/merged/all.genotypes.AF.0.05.vcf.gz'
	params:
		'imputation_files/posfiles/HRC.r1.GRCh37.autosomes.mac5.sites.AF.0.05.vcf.gz'
	shell:
		''' bcftools view -Oz -R {params} {input} > {output}
		bcftools index {output} '''
	
# generate ROC data for heterozygous SNP calls from imputation with SNP arrays as ground truth
rule accuracy:
	input:
		imputed = 'analysis/stitch/merged/all.genotypes.AF.0.05.vcf.gz',
		array_based = get_snp_array_vcf
	output:
		snps='analysis/roc/{ind}/{ind}.compare.snps.txt',
		roc='analysis/roc/{ind}/{ind}.roc.txt'
	params:
		ind_id = get_snp_array_id,
	shell:
		''' sh scripts/extract_snps.sh {input.imputed} {input.array_based} analysis/roc/{wildcards.ind} {params.ind_id} {wildcards.ind}
		Rscript scripts/accuracy.R {output.snps} > {output.roc} '''
		
# plot ROC for heterozygous SNP calls from imputation with SNP arrays as ground truth
rule plot_roc:
	input:
		get_roc_targets,
		outliers='tmp.outliers.txt' #tmp
#		outliers='analysis/ancestry/EA.PC1.PC2.outliers'
	output:
		plot='analysis/roc/plots/roc.pdf',
		stats='analysis/roc/plots/roc.stats'
	shell:
		'Rscript scripts/plot_roc.R {input} {output.plot} > {output.stats}'

#config['SNP_array_list']

rule all:
	input:
		get_targets


