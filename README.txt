This workflow uses STITCH to impute genotypes based on next gen sequencing data (in particular, ChIP-seq data). It allows multiple bam files from a given individual to be combined prior to imputation to increase imputation accuracy. Follow these steps to run it:

1. create the necessary conda environment from impute.env.yml. (ie, 'conda env create --file impute.env.yml')
2. update the config.yaml file if needed (this allows you to specify the reference genome and some parameters for STITCH). Usually nothing needs to be changed here.
3. add imputation reference files used by STICH and indicate their location in the config.yaml file. The files can be obtaine from: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html. The contents of this directory are as follows (in this case, the directory is named 'imputation_files':

imputation_files/1000GP_Phase3.sample

imputation_files/hapmap:
1000GP_Phase3_chr10.hap.gz  1000GP_Phase3_chr16.hap.gz  1000GP_Phase3_chr21.hap.gz  1000GP_Phase3_chr6.hap.gz
1000GP_Phase3_chr11.hap.gz  1000GP_Phase3_chr17.hap.gz  1000GP_Phase3_chr22.hap.gz  1000GP_Phase3_chr7.hap.gz
1000GP_Phase3_chr12.hap.gz  1000GP_Phase3_chr18.hap.gz  1000GP_Phase3_chr2.hap.gz   1000GP_Phase3_chr8.hap.gz
1000GP_Phase3_chr13.hap.gz  1000GP_Phase3_chr19.hap.gz  1000GP_Phase3_chr3.hap.gz   1000GP_Phase3_chr9.hap.gz
1000GP_Phase3_chr14.hap.gz  1000GP_Phase3_chr1.hap.gz   1000GP_Phase3_chr4.hap.gz
1000GP_Phase3_chr15.hap.gz  1000GP_Phase3_chr20.hap.gz  1000GP_Phase3_chr5.hap.gz

imputation_files/legend:
1000GP_Phase3_chr10.legend.gz  1000GP_Phase3_chr18.legend.gz  1000GP_Phase3_chr4.legend.gz
1000GP_Phase3_chr11.legend.gz  1000GP_Phase3_chr19.legend.gz  1000GP_Phase3_chr5.legend.gz
1000GP_Phase3_chr12.legend.gz  1000GP_Phase3_chr1.legend.gz   1000GP_Phase3_chr6.legend.gz
1000GP_Phase3_chr13.legend.gz  1000GP_Phase3_chr20.legend.gz  1000GP_Phase3_chr7.legend.gz
1000GP_Phase3_chr14.legend.gz  1000GP_Phase3_chr21.legend.gz  1000GP_Phase3_chr8.legend.gz
1000GP_Phase3_chr15.legend.gz  1000GP_Phase3_chr22.legend.gz  1000GP_Phase3_chr9.legend.gz
1000GP_Phase3_chr16.legend.gz  1000GP_Phase3_chr2.legend.gz
1000GP_Phase3_chr17.legend.gz  1000GP_Phase3_chr3.legend.gz

imputation_files/posfiles:
HRC.r1-1.GRCh37.wgs.mac5.sites.10.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.21.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.11.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.22.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.12.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.2.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.13.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.3.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.14.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.4.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.15.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.5.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.16.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.6.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.17.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.7.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.18.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.8.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.19.EUR_AF01.posfile  HRC.r1-1.GRCh37.wgs.mac5.sites.9.EUR_AF01.posfile
HRC.r1-1.GRCh37.wgs.mac5.sites.1.EUR_AF01.posfile   HRC.r1.GRCh37.autosomes.mac5.sites.AF.0.05.vcf.gz
HRC.r1-1.GRCh37.wgs.mac5.sites.20.EUR_AF01.posfile

4. update the metasheet.csv file. This is file contains comma-separated columns with the bam location in the first column and the individual the bam corresponds to in the second column. For example:

bam,individual
chips/align/InputV16/InputV16_unique.sorted.dedup.bam,VCaP16
chips/align/InputVca/InputVca_unique.sorted.dedup.bam,VCaP
chips/align/R1V16K27/R1V16K27_unique.sorted.dedup.bam,VCaP16
chips/align/R1VCK27/R1VCK27_unique.sorted.dedup.bam,VCaP

5. activate the conda environment ('conda activate impute')
6. you can submit the workflow to slurm using submit.sh. Or, if not using slurm, run the workflow as, for example, 'snakemake -s impute.snakefile -pr all 
