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

g <- ggplot2::ggplot(AVR.gt5.df,
                ggplot2::aes(x = Year_Collected,
                             group = AVR,
                             y = Detected)) +
  ggplot2::geom_point(pch = "X", # dodge position vertically
                      size = 3,
                      position = ggplot2::position_jitter(width = 1, height = 0.05)) +
  ggplot2::scale_y_continuous() +
  ggplot2::facet_wrap(~AVR,
                      nrow = 2) +
  ggplot2::geom_line(data = AVR.split,
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split,
                       ggplot2::aes(x = Year_Collected,
                                    ymin = pred - SE,
                                    ymax = pred + SE),
                       alpha = 0.2,
                       inherit.aes = FALSE,
                       color = "black") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 legend.position.inside = c(0.95, 0.85),
                 panel.spacing.x = ggplot2::unit(1, "lines")) +
  ggplot2::ylab("Predicted Prevalence of AVR Gene") +
  ggplot2::xlab("Year Sample Collected")

ggplot2::ggsave(g,
                filename = "../figures/AVR-SRR_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")

