#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem=46GB            # maximum memory per node
#SBATCH -p atlas            # standard node(s)
#SBATCH -J SRR_dump           # job name
#SBATCH --account="${ACCOUNT}" # account name
#SBATCH --array=1-18%48

# 329

# load modules
module load apptainer
export SINGULARITY_CACHEDIR=${TMPDIR}
export SINGULARITY_TMPDIR=${TMPDIR}

DOWNLOADDIR=/project/90daydata/${ACCOUNT}/Yulin

cd ${DOWNLOADDIR}

SRR=$(cat SRR_redo.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "${SRR}"

cd fastq_se

if [ ! -f ${SRR}_1.fastq.gz ]; then
	singularity exec ~/tools/sratools_latest.sif fasterq-dump ../SRR/${SRR} \
		--threads 6 \
		--temp ${TMPDIR} \
		--split-3 --skip-technical \
		--mem 4G
	# rm the singletons
	${SRR}.fastq

	# NOTE! bgzip does NOT error out if the input file cannot be found
	# this is fortunate, because this allows us to deal with the case where the are no paired reads
	bgzip -@6 --index ${SRR}.fastq

	# rm -rf ../SRR/${SRR}
fi
