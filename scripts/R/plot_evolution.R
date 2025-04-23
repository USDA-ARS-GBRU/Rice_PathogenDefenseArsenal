here::i_am("scripts/R/plot_evolution.R")

Country = NULL
meta <- read.table(here::here("inputs", "meta_combined.txt"), header = TRUE,
                   comment.char = "", sep = "\t") |>
  dplyr::mutate(Domestic = factor(Country == "USA", levels = c("TRUE", "FALSE")))

gene <- "Pita1"
gene <- "Piks"

seq <- NULL
labs <- read.table(here::here("sequences", "AVR_Trees", gene, "seqnames.txt"),
                   header = FALSE, col.names = c("NewName", "OldName"),
                   comment.char = "")

tr <- ape::read.tree(here::here("sequences", "AVR_Trees", gene, "pep.phy_phyml_tree.txt"))


if (gene == "Pita1") {
  tr <- ape::drop.tip(tr, c("seq00003", "seq00570", "seq00557"))
}

if (gene == "Piks") {
  tr <- ape::drop.tip(tr, c("seq00491"))
}

# get the unique sequences
tokeep <- read.table(here::here("sequences", "AVR_Trees", gene, "pep.aln.txt"),
                     header = FALSE, col.names = c("sample", "seq")) |>
  dplyr::distinct(seq, .keep_all = TRUE)

tr <- ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% tokeep$sample])

tr$tip.label <- labs$NewName[match(tr$tip.label, labs$OldName)]
cols <- c("red", "black")[meta$Domestic[match(tr$tip.label, meta$Sample)]]
cols[is.na(cols)] <- "grey"

tr$edge.length[tr$edge.length <= 0.000001] <- 0.0
plot(tr, type = "fan",
     # use.edge.length = FALSE,
     show.tip.label = FALSE
     )
ape::tiplabels(pch = 19, col = cols, cex = 0.5)

library(extrafont)

matchar_to_matint <- function(mat) {
  return(mat |>
           data.frame() |>
           dplyr::mutate(dplyr::across(dplyr::everything(),
                                       forcats::fct_infreq)) |>
           dplyr::mutate(dplyr::across(dplyr::everything(),
                                       as.integer)) |>
           data.frame() |>
           as.matrix())
}

SRR.meta <- read.table(here::here("inputs", "SRR_meta.txt"),
                       sep = "\t", header = TRUE, check.names = FALSE)
Wang.meta <- readRDS(here::here("inputs", "pathogen.RDS"))[, c(1:4)]

pep.Piks <- read.table(here::here("sequences",
                                  "AVR_Alignments",
                                  "Piks", "pep.aln.txt"),
                       header = FALSE, sep = "\t")

pep.Piks.mat <- do.call(rbind, pep.Piks$V2 |> strsplit(""))
rownames(pep.Piks.mat) <- pep.Piks$V1


pep.Pita <- read.table(here::here("sequences",
                                  "AVR_Alignments",
                                  "Pita1", "pep.aln.txt"),
                       header = FALSE, sep = "\t", comment.char = "")

pep.Pita.mat <- do.call(rbind, pep.Pita$V2 |> strsplit(""))
rownames(pep.Pita.mat) <- pep.Pita$V1

AA.Piks <- pep.Piks.mat[, c(46, 47, 48, 78)]
colnames(AA.Piks) <- paste0("codon_", c(46, 47, 48, 78))

AA.Pita <- pep.Pita.mat[, c(83, 119, 192, 207)]
colnames(AA.Pita) <- paste0("codon_",  c(83, 119, 192, 207))

AA.Piks <- pep.Piks.mat
AA.Pita <- pep.Pita.mat

# keep only sites with at least two amino acids
AA.Piks <- AA.Piks[, apply(AA.Piks, 2, \(x) length(table(x)) >= 2)]
AA.Pita <- AA.Pita[, apply(AA.Pita, 2, \(x) length(table(x)) >= 2)]

# convert into integer for clustering
AA.Pita.int <- matchar_to_matint(AA.Pita)
AA.Piks.int <- matchar_to_matint(AA.Piks)

# custom distance function (pointer)
catFuncPtr <- RcppXPtrUtils::cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) {
  return arma::accu(A != B)/static_cast<double>(A.n_cols);
                                     }",
                                     depends = c("RcppArmadillo"))


AA.Pita.dist <- parallelDist::parDist(AA.Pita.int, method="custom", func = catFuncPtr)
AA.Piks.dist <- parallelDist::parDist(AA.Piks.int, method="custom", func = catFuncPtr)

AA.Pita.dist.toplot <- AA.Pita.dist |> as.dist()
AA.Pita.dist.toplot[AA.Pita.dist.toplot > 0.20] <- 0.20
pheatmap::pheatmap(AA.Pita.dist.toplot)

AA.Piks.dist.toplot <- AA.Piks.dist |> as.dist()
AA.Piks.dist.toplot[AA.Piks.dist.toplot > 0.20] <- 0.20
pheatmap::pheatmap(AA.Piks.dist.toplot)

hclust(AA.Pita.dist, method = "complete") |>
  plot()

hclust(AA.Piks.dist, method = "complete") |>
  plot()
