#!/bin/bash

cd sequences/AVR_Trees/

GENE=Piks

cd ${GENE}
awk -v OFS="\t" '{printf("%s\tseq%05d\n", $1, NR)}' ../../AVR_Alignments/${GENE}/pep.aln.txt  > seqnames.txt
paste <(cut -f2 seqnames.txt) <(cut -f2 ../../AVR_Alignments/${GENE}/pep.aln.txt) | seqkit tab2fx > pep.aln.fasta
seqret -sequence pep.aln.fasta -outseq "phylip::pep.phy"
mpirun -np 20 phyml-mpi --datatype aa --input pep.phy --leave_duplicates --r_seed 1234
cd ..

GENE=Pita1

cd ${GENE}
awk -v OFS="\t" '{printf("%s\tseq%05d\n", $1, NR)}' ../../AVR_Alignments/${GENE}/pep.aln.txt  > seqnames.txt
paste <(cut -f2 seqnames.txt) <(cut -f2 ../../AVR_Alignments/${GENE}/pep.aln.txt) | seqkit tab2fx > pep.aln.fasta
seqret -sequence pep.aln.fasta -outseq "phylip::pep.phy"
mpirun -np 20 phyml-mpi --datatype aa --input pep.phy --leave_duplicates --r_seed 1234
cd ..

