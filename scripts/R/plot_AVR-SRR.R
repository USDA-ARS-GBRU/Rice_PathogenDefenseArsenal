setwd("../../inputs")

library(extrafont)

AVR_Primer <- NULL
Present <- NULL
F <- NULL
R <- NULL

AVR <- read.table("../outputs/AVR_primer_counts_in_SRR.txt.gz", check.names = FALSE)
meta <- read.table("SRR_meta.txt", header = TRUE, sep = "\t", check.names = FALSE) |>
  dplyr::rename("Year_Collected" = "Year Collected")

AVR.gt5 <- AVR >= 5L

AVR.gt5.df <- AVR.gt5 |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("ID") |>
  tidyr::pivot_longer(-1, names_to = "AVR_Primer", values_to = "Present") |>
  tidyr::separate(AVR_Primer, c("AVR", "Primer"), sep = "_") |>
  tidyr::pivot_wider(names_from = "Primer", values_from = Present) |>
  dplyr::mutate(Detected = `F` & `R`) |>
  dplyr::mutate(Detected = Detected * 1L) |>
  dplyr::left_join(meta, by = c("ID" = "SRR Sample ID"))

# estimate some smooths from mgcv
AVR.split <- split(AVR.gt5.df, AVR.gt5.df$AVR) |>
  lapply(\(x) mgcv::gam(cbind(Detected, 1 - Detected) ~ 1 + s(Year_Collected, bs = "cs"),
                        data = x,
                        method = "REML",
                        family = "binomial")) |>
  lapply(predict,
         newdata = data.frame(Year_Collected = 1950:2020),
         se.fit = TRUE,
         type = "response") |>
  lapply(\(x) data.frame(Year_Collected = 1950:2020,
                         pred = x$fit |> unname(),
                         SE = x$se.fit |> unname())) |>
  purrr::list_rbind(names_to = "AVR")

AVR <- AVR.gt5.df
g <- ggplot2::ggplot(AVR,
                     ggplot2::aes(x = Year_Collected,
                                  group = AVR,
                                  y = Detected)) +
  ggplot2::geom_point(pch = 21, # dodge position vertically
                      size = 1,
                      alpha = 0.7,
                      color = "darkgrey",
                      mapping = ggplot2::aes(fill = as.factor(Detected)),
                      position = ggplot2::position_jitter(width = 1, height = 0.1)) +
  ggplot2::facet_wrap(~AVR,
                      nrow = 2) +
  ggplot2::geom_line(data = AVR.split,
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split,
                       ggplot2::aes(x = Year_Collected,
                                    ymin = pmax(pred - SE, 0),
                                    ymax = pmin(pred + SE, 1)),
                       alpha = 0.2,
                       inherit.aes = FALSE,
                       color = "black") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 # legend.position = "none",
                 legend.position = "bottom",
                 # legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 # legend.position.inside = c(0.95, 0.85),
                 panel.spacing.x = ggplot2::unit(1, "lines")) +
  ggplot2::scale_y_continuous(name = "Predicted Prevalence of AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25)) +
  ggplot2::xlab("Year Sample Collected") +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             # values = viridisLite::plasma(2),
                             values = c("purple", "lightgrey"),
                             labels = c("Absent", "Present")) +
  ggplot2::ggtitle("AVR Presence/Absence in Pathogen Samples, over Time (NCBI Samples)")
g
ggplot2::ggsave(g,
                filename = "../figures/AVR-SRR_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")

AVR.Global <- AVR
AVR.split.Global <- AVR.split

save(AVR.Global,
     AVR.split.Global,
     file = "../outputs/source_data_AVR_PA_over_time_Global.RData")

