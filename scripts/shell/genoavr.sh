#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem=46GB            # maximum memory per node
#SBATCH -p atlas            # standard node(s)
#SBATCH -J ssr_geno           # job name
#SBATCH --account="${ACCOUNT}" # account name
#SBATCH --array=1-1

module load apptainer

cd /project/90daydata/${ACCOUNT}/Yulin/avr2

rm -f ${SSR}.txt

for SAMPLE in $(bcftools head ../pixy_new/combined.vcf.gz | fgrep -v "##" | cut -f10-);
	do
	echo ${SAMPLE}
	samtools view --with-header --bam ../alignment*/${SAMPLE}.bam ${REGION} | samtools addreplacerg -w -r ID:${SAMPLE} -r LB:${SAMPLE} -r SM:${SAMPLE} -o ${SAMPLE}_${SSR}.bam -
	samtools index ${SAMPLE}_${SSR}.bam
	# try agian with HC
	singularity exec ~/tools/gatk_4.6.1.0.sif gatk --java-options "-Xms7G -Xmx7G" HaplotypeCaller \
		--tmp-dir ${TMPDIR} \
		--native-pair-hmm-threads 1 \
		--reference ../Magor1_AssemblyScaffolds.fasta.gz \
		--input ${SAMPLE}_${SSR}.bam \
		--output ${SSR}_${SAMPLE}.vcf.gz \
		--bam-output ${SSR}_${SAMPLE}.bam \
		--bam-writer-type CALLED_HAPLOTYPES \
		--sample-ploidy 1 \
		--pileup-detection \
		--emit-ref-confidence BP_RESOLUTION \
		--indel-size-to-eliminate-in-ref-model 100 \
		--max-alternate-alleles 4 \
		--base-quality-score-threshold 15 \
		--min-base-quality-score 10 \
		--minimum-mapping-quality 10 \
		--mapping-quality-threshold-for-genotyping 10 \
		--pcr-indel-model NONE \
		--annotation-group StandardAnnotation \
		--annotation-group AS_StandardAnnotation \
		--annotation-group StandardHCAnnotation \
		--intervals ${REGION}
	# --include-non-variant-sites \
	singularity exec ~/tools/gatk_4.6.1.0.sif gatk --java-options "-Xms30G -Xmx36G" GenotypeGVCFs \
		--tmp-dir ${TMPDIR} \
		--variant ${SSR}_${SAMPLE}.vcf.gz \
		--reference ../Magor1_AssemblyScaffolds.fasta.gz \
		--intervals ${REGION}\
		--only-output-calls-starting-in-intervals \
		--annotation-group StandardAnnotation \
		--annotation-group AS_StandardAnnotation \
		--max-alternate-alleles 2 \
		--sample-ploidy 1 \
		--output ${SSR}_${SAMPLE}.call.vcf.gz
	echo -e ${SAMPLE} "\t" $(samtools faidx ../Magor1_AssemblyScaffolds.fasta.gz ${REGION} | \
		bcftools consensus --regions-overlap 0 --haplotype 1 --prefix ${SAMPLE} --samples ${SAMPLE} ${SSR}_${SAMPLE}.call.vcf.gz |
                seqkit seq --reverse --complement | \
		cut -f2 ) >> ${SSR}.txt
        rm ${SAMPLE}_${SSR}.bam*
done

#                 seqkit fx2tab | \
