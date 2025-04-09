#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem=46GB            # maximum memory per node
#SBATCH -p atlas            # standard node(s)
#SBATCH -J avr_geno           # job name
#SBATCH --account="${ACCOUNT}" # account name
#SBATCH --array=1-329%48

DOWNLOADDIR=/project/90daydata/${ACCOUNT}/Yulin

cd ${DOWNLOADDIR}

SRR=$(cat SRR.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "${SRR}"

seqkit locate --max-len-to-show 1 --max-mismatch 0 --threads 5 -f primers.fasta ./trim_reads*/${SRR}_trim_*.fastq.gz |
        tail -n+2 |
        cut -f2 |
        sort |
        uniq -c |
        awk '{print $2,$1}' > ./avr/${SRR}.txt
