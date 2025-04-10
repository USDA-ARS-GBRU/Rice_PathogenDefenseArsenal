#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB            # maximum memory per node
#SBATCH -p service            # standard node(s)
#SBATCH -J SRR_down           # job name
#SBATCH --account="${ACCOUNT}" # account name
#SBATCH --array=1-18

# load modules
module load apptainer
export SINGULARITY_CACHEDIR=${TMPDIR}
export SINGULARITY_TMPDIR=${TMPDIR}

DOWNLOADDIR=/project/90daydata/${ACCOUNT}/Yulin

cd ${DOWNLOADDIR}

SRR=$(cat SRR_redo.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

cd SRR
if [ -f ${SRR}_2.fastq.gz ]; then
	echo "FILE ALREADY DOWNLOADED: ${SRR}";
	exit 0;
fi

if [ ! -f ${SRR}_2.fastq.gz ]; then
	sleep $((RANDOM % 60))
	singularity exec ~/tools/sratools_latest.sif prefetch "$SRR" --max-size u
fi

# retry 2
if [ ! -d ${SRR} ]; then
        sleep $((RANDOM % 60))
        singularity exec ~/tools/sratools_latest.sif prefetch "$SRR" --max-size u
fi

# retry 3
if [ ! -d ${SRR} ]; then
        sleep $((RANDOM % 60))
        singularity exec ~/tools/sratools_latest.sif prefetch "$SRR" --max-size u
fi

# retry 4
if [ ! -d ${SRR} ]; then
        sleep $((RANDOM % 60))
        singularity exec ~/tools/sratools_latest.sif prefetch "$SRR" --max-size u
fi

# retry 5
if [ ! -d ${SRR} ]; then
        sleep $((RANDOM % 60))
        singularity exec ~/tools/sratools_latest.sif prefetch "$SRR" --max-size u
fi
