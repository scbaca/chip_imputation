#!/bin/bash
#
#SBATCH --job-name=impute
#SBATCH --output=out.impute.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH -t 10:00:00

srun snakemake -s impute.snakefile -j 100 -pr --rerun-incomplete --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads}" all

#srun snakemake -s impute.snakefile -j 100 -pr --rerun-incomplete --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads}" analysis/genotype/plots/gt.heatmap.pdf

