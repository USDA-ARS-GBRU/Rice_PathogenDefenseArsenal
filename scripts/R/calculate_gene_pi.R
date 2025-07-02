
# set here:: ------------------------------------------------------------------------------------------------------

here::i_am("scripts/R/calculate_gene_pi.R")




# set ignores -----------------------------------------------------------------------------------------------------

Country <- NULL


# load Arial font -------------------------------------------------------------------------------------------------

library(extrafont)
extrafont::loadfonts(device="postscript")

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

width <- 48

if (gene == "Pita1") {
  # size 22
  meta$Group <- ifelse(meta$Year < 1995,
                       "Pre 1995",
                       ifelse(meta$Year <= 2005,
                              "1995 - 2005",
                              "Post 2005")) |>
    factor(levels = c("Pre 1995", "1995 - 2005", "Post 2005"))
} else if (gene == "Piks") {
  # size 27
  meta$Group <- ifelse(meta$Year < 1990,
                       "Pre 1990",
                       ifelse(meta$Year <= 2010,
                              "1990 - 2010",
                              "Post 2010")) |>
    factor(levels = c("Pre 1990", "1990 - 2010", "Post 2010"))
}

base <- seq(from = 1, to = nrow(aln) - width, by = 6)

popgen <- split(meta$Sample, interaction(ifelse(meta$Domestic, "U.S.", "Global"), meta$Group, sep = "_")) |>
  lapply(\(ZZ)
    lapply(base,
                 \(x) seq(from = x, to = x + width)) |>
      lapply(\(x) genomicpp::rcpp_parallel_tajimas_d(aln[x, colnames(aln) %in% ZZ])) |>
      lapply(as.data.frame) |>
      lapply(t) |>
      lapply(as.data.frame) |>
      purrr::list_rbind() |>
      dplyr::mutate("base" = base, .before = 1) |>
      tibble::remove_rownames()) |>
  purrr::list_rbind(names_to = "Domestic_Group") |>
  tidyr::separate(Domestic_Group, c("Domestic", "Group"), sep = "_")

if (gene == "Pita1") {
  popgen$Group <- factor(popgen$Group,
                         levels = c("Pre 1995", "1995 - 2005", "Post 2005"))
} else if (gene == "Piks") {
  popgen$Group <- factor(popgen$Group,
                         levels = c("Pre 1990", "1990 - 2010", "Post 2010"))
}

popgen.long <- popgen |>
  tidyr::pivot_longer(-c(1:3), names_to = "Statistic", values_to = "value")

ggplot2::ggplot(popgen.long |>
                  dplyr::filter(Statistic %in% c("pi", "D")),
                ggplot2::aes(x = base,
                             y = value,
                             color = Domestic)) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(Statistic~Group, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("Global" = "lightblue", "U.S." = "salmon2")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "bottom",
                 panel.spacing.x = ggplot2::unit(1, "lines"))

# try to do some bootstrapping....
bases <- split(meta$Sample, interaction(ifelse(meta$Domestic, "U.S.", "Global"), meta$Group, sep = "_"))
bases <- lapply(bases, \(x) Filter(\(y) y %in% colnames(aln), x))
minsz <- min(sapply(bases, length))
minsz

set.seed(1234)

# generate 1000 samples from each group with replacement
samples <- lapply(bases,
       \(x) lapply(1:1000,
              \(i) sample(x, size = minsz, replace = TRUE)))

# samples <- lapply(bases, list)

samples.piD <- lapply(samples,
                      \(x) lapply(x,
                                  \(seqs) {
                                    aln.temp <- t(aln[, seqs]) |>
                                      as.data.frame() |>
                                      lapply(as.factor) |>
                                      lapply(as.integer) |>
                                      do.call(what = cbind) |>
                                      t();
                                    aln.temp[] <- aln.temp - 1L;
                                    genomicpp::rcpp_parallel_tajimas_d(aln.temp);
                                    }) |>
                        lapply(\(x) x[c("pi", "D")]) |>
                        lapply(t) |>
                        do.call(what = "rbind"))

samples.piD.long <- lapply(samples.piD, \(x) apply(x, MARGIN = 2, quantile, c(0.025, 0.5, 0.975))) |>
  lapply(as.data.frame) |>
  lapply(tibble::rownames_to_column, "perc") |>
  lapply(tidyr::pivot_longer, 2:3, names_to = "Statistic", values_to = "value") |>
  purrr::list_rbind(names_to = "Domestic_Group") |>
  tidyr::separate(Domestic_Group, c("Domestic", "Group"), sep = "_") |>
  tidyr::pivot_wider(names_from = "perc", values_from = "value")

if (gene == "Pita1") {
  samples.piD.long$Group <- factor(samples.piD.long$Group,
                         levels = c("Pre 1995", "1995 - 2005", "Post 2005"))
} else if (gene == "Piks") {
  samples.piD.long$Group <- factor(samples.piD.long$Group,
                         levels = c("Pre 1990", "1990 - 2010", "Post 2010"))
}

samples.piD.long

# plotting code... probably not useful

# ggplot2::ggplot(samples.piD.long,
#                 ggplot2::aes(x = Group,
#                              y = `50%`,
#                              color = Domestic,
#                              group = Domestic)) +
#   ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.2)) +
#   ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.2)) +
#   ggplot2::facet_wrap(~Statistic, nrow = 2, scales = "free_y") +
#   ggplot2::scale_color_manual(values = c("Global" = "lightblue", "U.S." = "salmon2")) +
#   ggplot2::theme_minimal() +
#   ggplot2::geom_errorbar(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
#                        width = 0.2,
#                        position = ggplot2::position_dodge(width = 0.2)) +
#   ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
#                  legend.position = "bottom",
#                  panel.spacing.x = ggplot2::unit(1, "lines")) +
#   ggh4x::facetted_pos_scales(y = list(
#     Statistic == "D" ~ ggplot2::scale_y_continuous(limits = c(-4, 0)),
#     Statistic == "pi" ~  ggplot2::scale_y_continuous(limits = c(0, 0.02))
#     )
#   ) +
#   ggplot2::ylab("Statistic Value")
