here::i_am("scripts/R/plot_evolution.R")

Country = NULL
meta <- read.table(here::here("inputs", "meta_combined.txt"), header = TRUE,
                   comment.char = "", sep = "\t") |>
  dplyr::mutate(Domestic = ifelse(Country == "USA", TRUE, FALSE))

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

plot(tr, use.edge.length = FALSE)

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

load(here::here("outputs",
                "source_data_AVR_PA_over_time_Global.RData"))

library(ape)
####
# read in a coding sequence MSA
plot_msa_network <- function(fa = here::here("sequences", "AVR_Alignments", "Pita1", "cds.aln.fasta"),
                             img = here::here("figures", "Pita1.png"),
                             metadata = meta) {
  gene <- fs::path_ext_remove(basename(img)) |>
    stringr::str_replace_all("haplotype\\_network\\_",
                             "")

  seqs <- adegenet::fasta2DNAbin(fa)

  # drop individuals that failed the "PCR" test
  bad <- subset(AVR.Global,
                AVR == paste0("AVR-", gene) &
                  Detected == 0,
                select = ID,
                drop = TRUE)

  seqs <- seqs[!dimnames(seqs)[[1]] %in% bad, ]

  saveRDS(seqs, file = here::here("outputs", paste0(gene, "_dnabin.RDS")))

  # subset seqs to drop individuals that didn't actually have the gene detected...

  # extract all of the individual haplotype sequences
  haplo <- pegas::haplotype(seqs,
                            trailingGapsAsN = FALSE,
                            strict = FALSE)

  # rename haplotype sequences from roman numerals to VXXX
  attr(haplo, "dimnames")[[1]] <- paste0("V",
                                         stringr::str_pad(seq_along(attr(haplo, "dimnames")[[1]]),
                                                          width = 3, pad = "0"))

  # get the groupings
  haplo.index.list <- attr(haplo, "index")
  names(haplo.index.list) <- dimnames(haplo)[[1]]


  # make the group assignments for putting in to a supp table
  supp <- data.frame("Sample" = dimnames(seqs)[[1]]) |>
    dplyr::left_join(lapply(haplo.index.list,
                            \(x) data.frame("Sample" = dimnames(seqs)[[1]][x])) |>
                       purrr::list_rbind(names_to = "hap"),
                     by = "Sample") |>
    dplyr::left_join(metadata[, c("Sample", "Domestic")],
                     by = "Sample") |>
    dplyr::mutate(Domestic = tidyr::replace_na(ifelse(Domestic, "Yes", "No"),
                                               "No")) |>
    dplyr::left_join(adegenet::DNAbin2genind(seqs) |>
                       adegenet::genind2df() |>
                       tibble::rownames_to_column("Sample"),
                     "Sample")

  if (gene == "Piks") {
    fname = here::here("tables_staged", "Dataset_S5.csv")
  } else {
    fname = here::here("tables_staged", "Dataset_S4.csv")
  }

  write.table(supp,
              file = fname,
              row.names = FALSE,
              col.names = TRUE,
              sep = ",",
              quote = FALSE,
              na = "-")

  # check if any of the samples from the US
  # RED = only US
  # BLUE = only non-US
  # BLACK = combination
  cols <- lapply(haplo.index.list, \(x) dimnames(seqs)[[1]][x]) |>
    sapply(\(x) x %in% metadata$Sample[metadata$Domestic]) |>
    sapply(\(x) ifelse(all(x), "salmon2",
                       ifelse(all(!x), "lightblue", "orchid")))

  # get the pairwise distances
  hnd <- ape::dist.dna(haplo, model = "N")

  # construct the haplotype network
  hnp <- pegas::haploNet(haplo, hnd, getProb = FALSE)

  # figure out labels
  # red if haplotype is found in US, false otherwise

  # reorder items to match the haplotype network
  cols <- cols[attr(hnp, "labels")]
  sizes <- attr(hnp, "freq")
  names(sizes) <- attr(hnp, "labels")

  p <- plot(hnp,
       size = sqrt(attr(hnp, "freq")),
       col = cols,
       show.mutation=1, scale.ratio=1.2,cex=0.6,threshold=0,lwd=0.8, fast = TRUE)

  # p

  # modified version of pegas:::as.igraph.haploNet
  # because it wasn't working correctly otherwise
  as.igraph.haploNet <- function(x, directed = FALSE, use.labels = TRUE,
           altlinks = TRUE, ...)
  {
    y <- x[, 1:2]
    if (altlinks) y <- rbind(y, attr(x, "alter.links")[, 1:2])
    y <-
      if (use.labels) matrix(attr(x, "labels")[y], ncol = 2)
    else y - 1L
    igraph::graph_from_edgelist(y, directed = directed, ...)
  }

  # convert the haplotype network to an igraph object
  ig <- as.igraph.haploNet(hnp, directed = FALSE, altlinks = FALSE)

  # get the names of the vertices in the igraph
  nm <- igraph::vertex.attributes(ig)$name

  # set the sizes of the nodes to be a minimum of 11
  # ie, number of times the sequence is seen plus 10
  # and then rescale by 5x (this seemed to work well for this dataset
  igraph::V(ig)$size <- (sizes[nm] + 10)*6
  igraph::V(ig)$size2 <- (sizes[nm] + 10)*6

  # set the seed for the (rotational aspect? of the) layout
  set.seed(1234)
  lay <- igraph::layout_with_fr(ig)

  cols <- lapply(haplo.index.list, \(x) dimnames(seqs)[[1]][x]) |>
    lapply(\(x) c(sum(x %in% metadata$Sample[metadata$Domestic]),
                  sum(!x %in% metadata$Sample[metadata$Domestic])))

  cols <- cols[attr(hnp, "labels")]

  # rescale layout to have equal height and width
  lay <- scale(lay)*15

  # salmon2 is US
  # lightblue is World

  png(img, width = 7, height = 7, units = "in", res = 600, pointsize = 7)
  par(bg = NA)
  par(mar = c(0, 0, 0, 0))
  plot(ig,
       xlim = range(lay[, 1]),
       ylim = range(lay[, 2]),
       layout = lay,
       vertex.shape = "pie",
       vertex.pie = cols[nm],
       vertex.pie.color = lapply(cols[nm], \(x) c("salmon2", "lightblue")[x > 0]),
       # vertex.pie.density = 1,
       rescale = FALSE,
       # asp = 0,
       vertex.label = NA, family = "Arial")

  plot(ig,
       xlim = range(lay[, 1]),
       ylim = range(lay[, 2]),
       layout = lay,
       vertex.size = 0,
       rescale = FALSE,
       vertex.label.cex = 0.85,
       vertex.label.color = "black",
       vertex.label.dist = 12,
       vertex.label.degree = -pi/2,
       add = TRUE, family = "Arial")


  if (gene == "Piks") {
    legend("bottomright",
           legend = c("U.S.", "Global"),
           col = c("salmon2", "lightblue"),
           pch = 19, inset = c(0.05, 0.05),
           cex = 3)
  }

  dev.off()
  }

plot_msa_network(fa = here::here("sequences", "AVR_Alignments", "Pita1", "cds.aln.fasta"),
                 img = here::here("figures", "haplotype_network_Pita1.png"))

plot_msa_network(fa = here::here("sequences", "AVR_Alignments", "Piks", "cds.aln.fasta"),
                 img = here::here("figures", "haplotype_network_Piks.png"))

# https://stackoverflow.com/questions/25360248/arrange-multiple-32-png-files-in-a-grid
# plot1 <- png::readPNG(here::here("figures", "haplotype_network_Pita1.png"))
# plot2 <- png::readPNG(here::here("figures", "haplotype_network_Piks.png"))
#
# tmp <- gridExtra::arrangeGrob(grid::rasterGrob(plot1), grid::rasterGrob(plot2), nrow = 1)
#
# ggplot2::ggsave(here::here("figures", "Figure_4.png"), tmp, width = 14, height = 7)
#
