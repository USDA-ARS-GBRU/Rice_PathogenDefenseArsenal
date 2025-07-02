here::i_am("scripts/R/plot_gpixy.R")
library(extrafont)


NAME <- NULL
START <- NULL
type <- NULL

regions <- read.table(here::here("outputs", "interval_popgen.txt.gz"),
                      header = TRUE, sep = "\t")

# look at our genes:
regions |>
  dplyr::filter(grepl("_4881_", NAME) |
                  grepl("_5074_", NAME) |
                  grepl("_10623_", NAME) |
                  grepl("_11753_", NAME) |
                  grepl("_12522_", NAME) |
                  grepl("_12603_", NAME))

quantile(regions$D, 0:20/20)

# ecdf(log10(ifelse(regions$pi == 0,
#                   min(regions$pi[regions$pi > 0])/10,
#                   regions$pi))) |>
#   plot(xlim = c(-8, 0), verticals = TRUE)
#
#
# ecdf(log10(ifelse(regions$pi == 0,
#                   min(regions$pi[regions$pi > 0])/10,
#                   regions$pi))) |>
#   plot(xlim = c(-8, 0),
#        verticals = TRUE)

regions <- regions |>
  tidyr::separate(NAME, into = c("type", "extra"), extra = "merge",
                  remove = FALSE)

ecdf.D <- subset(regions, D < 5 & STOP - START + 1L > 50L & type == "exon", select = D, drop = TRUE) |>
  sort()

#
# ecdfs.D <- split(regions$D, regions$type)
#
# D_obs <- clipr::read_clip() |> as.double()
# sapply(D_obs, \(x) mean(x < ecdfs.D$exon)) |>
#   clipr::write_clip()
#
# ecdfs.pi <- split(regions$pi, regions$type)
#
# pi_obs <- clipr::read_clip() |> as.double()
# sapply(pi_obs, \(x) mean(x > ecdfs.pi$exon)) |>
#   clipr::write_clip()

g1 <- ggplot2::ggplot(data = subset(regions, D < 5 & STOP - START + 1L > 50L),
                ggplot2::aes(x = START/1e6,
                             y = D,
                             fill = type)) +
  ggplot2::geom_point(alpha = 0.5, size = 2, color = "black", shape = 21) +
  ggplot2::geom_smooth(ggplot2::aes(color = type,
                                    linetype = type)) +
  ggplot2::facet_grid(~CHR, space = "free_x", scales = "free_x",
                      labeller = ggplot2::labeller(CHR = function(x) {
                        x <- as.character(x);
                        stringr::str_replace_all(x, "\\_", " ") |>
                          stringr::str_replace_all("Chromosome", "Chr.");
                      })) +
  ggplot2::theme_minimal() +
  ggplot2::xlab("Position on Chromosome (Mb)") +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 # legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.5),
                                                           color = "black"),
                 # legend.position.inside = c(0.94, 0.88),
                 panel.spacing.x = ggplot2::unit(0.5, "lines")) +
  ggplot2::ylab(latex2exp::TeX("Tajima's $\\textit{D}$"))

# ggplot2::ggsave(g1,
#                 filename = "../figures/Tajimas-D_over_time.png",
#                 width = 6,
#                 height = 2,
#                 dpi = 600,
#                 scale = 3,
#                 units = "in")

g2 <- ggplot2::ggplot(data = subset(regions, D < 5 & STOP - START + 1L > 50L),
                ggplot2::aes(x = START/1e6,
                             y = pi,
                             fill = type)) +
  ggplot2::geom_point(alpha = 0.5, size = 2, color = "black", shape = 21) +
  ggplot2::geom_smooth(ggplot2::aes(color = type,
                                    linetype = type)) +
  ggplot2::facet_grid(~CHR, space = "free_x", scales = "free_x") +
  ggplot2::theme_minimal() +
  ggplot2::geom_abline(intercept = mean(regions$D),
                       slope = 0,
                       col = "red", linetype = "dashed") +
  ggplot2::scale_y_continuous(name = latex2exp::TeX("$\\hat{\\pi}$ (Nucleotide Diversity)"),
                              trans = scales::pseudo_log_trans(base = 10,
                                                             sigma = 0.0001),
                              limits = c(0, 0.15),
                              breaks = c(0, 0.0001, 0.0005, 0.001,
                                         0.005, 0.01, 0.02, 0.03, 0.04,
                                         0.06, 0.08, 0.12),
                              labels = scales::label_number(drop0trailing=TRUE)) +
  ggplot2::xlab("Position on Chromosome (Mb)") +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 # legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.3),
                                                           color = "black"),
                 # legend.position.inside = c(0.94, 0.88),
                 panel.spacing.x = ggplot2::unit(0.5, "lines"),
                 strip.text = ggplot2::element_blank())

# ggplot2::ggsave(g2,
#                 filename = "../figures/pi_over_time.png",
#                 width = 6,
#                 height = 2,
#                 dpi = 600,
#                 scale = 3,
#                 units = "in")

f5.AB <- ggpubr::ggarrange(g1 + ggpubr::rremove("xlab") + ggpubr::rremove("x.ticks") + ggpubr::rremove("x.text"), g2,
                           common.legend = TRUE, legend = "none",ncol = 1,
                           heights = c(1, 1),
                           align = "v",
                           labels = c("A", "B"),
                           font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

g3 <- ggplot2::ggplot(data = subset(regions, D < 5 & STOP - START + 1L > 50L),
                ggplot2::aes(x =  pi,
                             color = type,
                             group = type,
                             linetype = type)) +
  ggplot2::stat_ecdf(geom = "step") +
  ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10,
                                                               sigma = 0.0001),
                              limits = c(0, 0.15),
                              breaks = c(0, 0.0001, 0.0005, 0.001,
                                         0.005, 0.01, 0.02, 0.03, 0.04,
                                         0.06, 0.08, 0.12),
                              labels = scales::label_number(drop0trailing = TRUE),
                              name = latex2exp::TeX("$\\hat{\\pi}$ (Nucleotide Diversity)")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "none", # "inside",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.3),
                                                           color = "black"),
                 # legend.position.inside =  c(0.80, 0.2),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 panel.spacing.x = ggplot2::unit(1, "lines"),
                 legend.title = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(name = latex2exp::TeX("\\overset{Cumulative Number of Features with}{$\\hat{\\pi}$ (Nucleotide Diversity) Less Than $\\textit{x}$}"),
                              position = "right")
# ggplot2::ggsave(g3,
#                 filename = "../figures/pi_ecdf.png",
#                 width = 5,
#                 height = 3,
#                 dpi = 600,
#                 scale = 2.2,
#                 units = "in")

f5.ABC <- ggpubr::ggarrange(f5.AB, g3,
                           common.legend = TRUE, legend = "top",
                           legend.grob = ggpubr::get_legend(g3),
                           ncol = 2,
                           widths = c(1.66, 1),
                           labels = c("", "C"),
                           font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

ggplot2::ggsave(plot = f5.ABC, filename = "../figures/Figure_5.png",
                dpi = 600, width = 7.5, height = 2.5, scale = 3)

