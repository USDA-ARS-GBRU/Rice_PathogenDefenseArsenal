
# set here:: ------------------------------------------------------------------------------------------------------

here::i_am("scripts/R/calculate_gene_pi.R")



# set ignores -----------------------------------------------------------------------------------------------------

Country <- NULL

# get sample metadata ---------------------------------------------------------------------------------------------

meta <- read.table(here::here("inputs", "meta_combined.txt"), header = TRUE,
                   comment.char = "", sep = "\t") |>
  dplyr::mutate(Domestic = ifelse(Country == "USA", TRUE, FALSE))

# set gene names --------------------------------------------------------------------------------------------------
gene <- "Pita1"
gene <- "Piks"

ok <- dimnames(readRDS(file = here::here("outputs", paste0(gene, "_dnabin.RDS"))))[[1]]

aln.raw <- read.table(here::here("sequences",
                                 "AVR_Alignments",
                                 gene, "cds.aln.txt"),
                      header = FALSE, sep = "\t", comment.char = "")

aln.mat <- do.call(rbind, aln.raw$V2 |>
                     strsplit(""))
aln.mat[aln.mat == "-"] <- NA_character_

rownames(aln.mat) <- aln.raw$V1

aln.mat <- t(aln.mat)
aln.mat <- aln.mat[, colnames(aln.mat) %in% ok]

# set missing val expected by genomicpp
aln <- matrix(-9L, nrow = nrow(aln.mat), ncol = ncol(aln.mat),
              dimnames = dimnames(aln.mat))

aln.new <- t(aln.mat) |>
  as.data.frame() |>
  lapply(as.factor) |>
  lapply(as.integer) |>
  do.call(what = cbind) |>
  t()

aln.new[] <- aln.new - 1L
aln[] <- aln.new

width <- 50

split(meta$Sample, meta$Domestic) |>
  lapply(\(ZZ) {
    pl <- lapply(seq(from = 1, to = nrow(aln) - width, by = 5),
           \(x) seq(from = x, to = x + width)) |>
      lapply(\(x) genomicpp::rcpp_parallel_tajimas_d(aln[x, colnames(aln) %in% ZZ])) |>
      lapply(as.data.frame) |>
      lapply(t) |>
      lapply(as.data.frame) |>
      purrr::list_rbind()

    plot(pl$pi, type = "l")
    plot(pl$D, type = "l")
})

