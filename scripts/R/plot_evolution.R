setwd("../../inputs")

library(extrafont)

SRR.meta <- read.table("SRR_meta.txt", sep = "\t", header = TRUE, check.names = FALSE)
Wang.meta <- readRDS("pathogen.RDS")[, c(1:4)]

pep.Piks <- read.table("../sequences/AVR_Alignments/Piks/pep.aln.txt",
                       header = FALSE, sep = "\t")

pep.Piks.mat <- do.call(rbind, pep.Piks$V2 |> strsplit(""))
rownames(pep.Piks.mat) <- pep.Piks$V1


pep.Pita <- read.table("../sequences/AVR_Alignments/Pita1/pep.aln.txt",
                       header = FALSE, sep = "\t", comment.char = "")

pep.Pita.mat <- do.call(rbind, pep.Pita$V2 |> strsplit(""))
rownames(pep.Pita.mat) <- pep.Pita$V1

AA.Piks <- pep.Piks.mat[, c(46, 47, 48, 78)]
colnames(AA.Piks) <- paste0("codon_", c(46, 47, 48, 78))

AA.Pita <- pep.Pita.mat[, c(83, 119, 192, 207)]
colnames(AA.Pita) <- paste0("codon_",  c(83, 119, 192, 207))


