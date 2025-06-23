here::i_am("scripts/R/plot_evolution.R")

library(extrafont)
# extrafont::font_import(prompt = FALSE)
extrafont::loadfonts(device="postscript")


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

# load in sequences with known missing status

load(here::here("outputs",
                "source_data_AVR_PA_over_time_Global.RData"))

bad <- subset(AVR.Global,
              AVR == paste0("AVR-", gene) &
                Detected == 0,
              select = ID,
              drop = TRUE)

fa <-  here::here("sequences", "AVR_Alignments", gene, "pep.aln.fasta")

seqs <- seqinr::read.alignment(file = fa, format = "fasta", seqtype = "AA")

# drop individuals that failed the "PCR" test
bad <- subset(AVR.Global,
              AVR == paste0("AVR-", gene) &
                Detected == 0,
              select = ID,
              drop = TRUE)

# figure out which sequences to keep
keep <- which(!seqs$nam %in% bad)

# drop the ones we dont want
# subset seqs to drop individuals that didn't actually have the gene detected...
seqs$seq <- seqs$seq[keep]
seqs$nam <- seqs$nam[keep]
seqs$nb <- length(keep)

seqs.mat <- as.matrix(seqs) |>
  toupper()


if (gene == "Pita1") {
  seqs.mat.crit <- seqs.mat[, c(83L, 119L, 192L, 207L)]
}

if (gene == "Piks") {
  seqs.mat.crit <- seqs.mat[, c(46L, 47L, 48L, 78L)]
}

# count the unique amino acid haplotypes
seqs.mat.crit |>
  apply(1, paste0, collapse = "") |> table() |> sort(decreasing = TRUE)

# combine with year information and plot over time???
# really not sure how this is gonna work out but can try I guess
seqs.df1 <- seqs.mat.crit |>
  apply(1, paste0, collapse = "") |>
  as.data.frame() |>
  setNames(c("AA")) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::left_join(meta[, c("Sample", "Year", "Domestic")], by = "Sample") |>
  dplyr::filter(!is.na(Year)) |>
  dplyr::arrange(Year) |>
  dplyr::mutate(AA = factor(AA, levels = unique(AA)))

g1 <- ggplot2::ggplot(seqs.df1,
                ggplot2::aes(x = Year,
                             y = Domestic,
                             color = Domestic,
                             group = Domestic)) +
  ggplot2::geom_point(pch = "|", size = 3) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~AA, ncol = 1, strip.position = "left") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_blank(),
                 legend.position = "inside",
                 legend.position.inside = c(0.1, 0.05),
                 legend.title = ggplot2::element_blank(),
                 strip.text = ggplot2::element_blank(),
                 panel.spacing.y = ggplot2::unit(0.2, "lines"),
                 legend.background = ggplot2::element_rect(fill = "white",
                                                           color = "black"),
                 text = ggplot2::element_text(size = 14, family = "Arial", color = "black")) +
  ggplot2::scale_color_manual(name = "",
                              values = c("lightblue", "salmon2"),
                              labels = c("Global", "U.S."),
                              guide = ggplot2::guide_legend(reverse = TRUE))

seqs.df2 <- seqs.mat.crit |>
  as.data.frame() |>
  dplyr::rename_with(.fn = \(x) paste0("codon_", sprintf("%03d", as.integer(x)))) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::left_join(seqs.df1, by = "Sample") |>
  dplyr::filter(!is.na(Year)) |>
  dplyr::select(2:6) |>
  dplyr::distinct() |>
  dplyr::arrange(as.integer(AA)) |>
  dplyr::mutate(dplyr::across(1:4, \(x) as.integer(factor(x, levels = unique(x))))) |>
  tidyr::pivot_longer(-5, names_to = "codon", values_to = "AA_val")

seqs.df2.right <- seqs.mat.crit |>
  as.data.frame() |>
  dplyr::rename_with(.fn = \(x) paste0("codon_", sprintf("%03d", as.integer(x)))) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::left_join(seqs.df1, by = "Sample") |>
  dplyr::filter(!is.na(Year)) |>
  dplyr::select(2:6) |>
  dplyr::distinct() |>
  dplyr::arrange(as.integer(AA)) |>
  tidyr::pivot_longer(-5, names_to = "codon", values_to = "AA_letter")

seqs.df2 <- seqs.df2 |>
  dplyr::left_join(seqs.df2.right, by = c("AA", "codon")) |>
  dplyr::mutate(AA_val = factor(AA_val))

g2 <- ggplot2::ggplot(seqs.df2,
                ggplot2::aes(x = codon,
                             y = 1,
                             fill = codon,
                             alpha = AA_val)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~AA, ncol = 1, strip.position = "right") +
  ggplot2::theme_minimal() +
  ggplot2::geom_text(ggplot2::aes(x = codon, y = 1, label = AA_letter),
                     alpha = 1,
                     size = 8,
                     # fontface = "bold",
                     family = "mono") +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "inside",
                 legend.position.inside = c(0.05, 0.05),
                 legend.title = ggplot2::element_blank(),
                 strip.text = ggplot2::element_blank(),
                 panel.spacing.y = ggplot2::unit(0.2, "lines"),
                 # panel.background = ggplot2::element_rect(fill = NA, color = "black"),
                 legend.background = ggplot2::element_rect(fill = "white",
                                                           color = "black"),
                 text = ggplot2::element_text(size = 14, family = "Arial", color = "black")) +
  ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(ifelse(gene == "Pita1", "Set1", "Set2"),
                                                               n = length(unique(seqs.df2$codon))),
                             guide = "none") +
  ggplot2::scale_alpha_manual(values = seq(0.8, 0.3, length.out = length(unique(seqs.df2$AA_val))),
                              guide = "none") +
  ggplot2::scale_x_discrete(labels = \(x) substr(x, 7, 9))

gp <- ggpubr::ggarrange(g1, g2, widths = c(3, 1), align = "h",
                        label.x = c(0, -0.1),
                        labels = c("A", "B"))

ggplot2::ggsave(gp,
                filename = here::here("figures", paste0("Figure_5_", gene, ".png")),
                width = 7,
                height = 7,
                dpi = 600,
                scale = 1,
                units = "in")

