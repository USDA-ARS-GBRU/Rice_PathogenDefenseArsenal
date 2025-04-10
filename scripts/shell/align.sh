#!/bin/bash
#SBATCH --cpus-per-task 48
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J bwt_algn               # job name
#SBATCH --array=1-311%32
#SBATCH --account="${ACCOUNT}" # account name

cd /project/90daydata/${ACCOUNT}/Yulin

mkdir -p alignments

cd ./alignments

PREFIX=$(cat ../SRR_pe.txt  | sed -n ${SLURM_ARRAY_TASK_ID}p)

if [[ -f ${PREFIX}.sam ]] ; then
    echo file exists.
    exit 0
fi

if [[ -f ${PREFIX}_MAPQ_pass ]] ; then
    echo file exists.
    exit 0
fi

echo "DEBUG: " $(date)

bwa-mem2 mem -t46 -o ${PREFIX}.sam ../Magor1_AssemblyScaffolds.fasta.gz ../trim_reads/${PREFIX}_trim_1.fastq.gz ../trim_reads/${PREFIX}_trim_2.fastq.gz

echo "DEBUG: " $(date)

samtools view --with-header -O sam --min-MQ 10 ${PREFIX}.sam | samtools sort -@ 48 -m 7G -O bam -l 9 -o ${PREFIX}.bam --write-index

echo "DEBUG: " $(date)

rm ${PREFIX}.sam
