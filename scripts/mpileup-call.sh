#!/bin/bash
#SBATCH -N 1                      # number of nodes
#SBATCH -n 48                     # n processor core(s) per node X 2 threads per core
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J bcf_call               # job name
#SBATCH --account="${ACCOUNT}" # account name

cd /project/90daydata/${ACCOUNT}/Yulin

mkdir -p calls

cd calls

bcftools mpileup \
    --output-type u \
    --fasta-ref ../Magor1_AssemblyScaffolds.fasta.gz \
    --max-depth 1000 \
    --ignore-RG \
    --threads 46 \
    --annotate FORMAT/AD \
    --full-BAQ \
    --min-MQ 10 \
    --min-BQ 15 \
    --min-ireads 5 \
    --regions ${CHR} \
    ../alignments/*.bam ../alignments_se/*.bam |
bcftools call \
    --multiallelic-caller \
    --ploidy 1 \
    --output-type b \
    --output calls.${CHR}.bcf \
    --threads 46 \
    --annotate FORMAT/GQ,FORMAT/GP \
    --write-index
