here::i_am("scripts/R/plot_evolution.R")


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

AA.Piks <- AA.Piks[, apply(AA.Piks, 2, \(x) length(table(x)) >= 2)]
AA.Pita <- AA.Pita[, apply(AA.Pita, 2, \(x) length(table(x)) >= 2)]

AA.Pita.int <- matchar_to_matint(AA.Pita)
AA.Piks.int <- matchar_to_matint(AA.Piks)

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

