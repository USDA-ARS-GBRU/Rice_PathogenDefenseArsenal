#!/bin/bash
#SBATCH --cpus-per-task 12                 # n processor core(s) per node X 2 threads per core
#SBATCH --mem=93GB            # maximum memory per node
#SBATCH -p atlas              # standard node(s)
#SBATCH -J fastp              # job name
#SBATCH --array=1-20%48
#SBATCH --account="${ACCOUNT}" # account name

cd /project/90daydata/${ACCOUNT}/Yulin

# make sure the final output folder is available
mkdir -p trim_reads_se
mkdir -p trim_reads_se_metrics

cd ./trim_reads_se

# skip if this one has already been processed

PREFIX=$(cat ../SRR_se.txt  | sed -n ${SLURM_ARRAY_TASK_ID}p)

if [[ -f ${PREFIX}_trim_1.fastq.gz ]] ; then
    echo file exists.
    exit 0
fi

# now run the read trimmer
/project/${ACCOUNT}/grant/daniel/tools/fastp \
    --in1 ../fastq_se/${PREFIX}.fastq.gz \
    --out1 ${PREFIX}_trim_1.fastq.gz \
    --thread 12 \
    --trim_front1 1 \
    --trim_tail1 0 \
    --cut_front \
    --cut_front_window_size 3 \
    --cut_front_mean_quality 15 \
    --cut_tail \
    --cut_tail_window_size 3 \
    --cut_tail_mean_quality 15 \
    --n_base_limit 5 \
    --average_qual 10 \
    --length_required 45 \
    --dedup \
    --dup_calc_accuracy 6 \
    --reads_to_process 450000000 \
    -h ../trim_reads_se_metrics/${PREFIX}.fastp.html \
    -j ../trim_reads_se_metrics/${PREFIX}.fastp.json
