here::i_am("scripts/R/plot_figureS1.R")

AVR.Global <- AVR.split.Global <- AVR.split.USA <- AVR.USA <- NULL
Rgenes <- Rgenes.loess <-  NULL
Rgenes_AVR_over_time <- Rgenes_AVR_over_time2 <- NULL
Rgenes_AVR_over_time.lag <- Rgenes_AVR_over_time.lag2 <- NULL

AVR <- R <- proportion <- NULL

load(file = here::here("outputs", "source_data_Rgenes_AVR_over_time.RData"))
load(file = here::here("outputs", "source_data_AVR_PA_over_time_USA.RData"))
load(file = here::here("outputs", "source_data_AVR_PA_over_time_Global.RData"))
load(file = here::here("outputs", "source_data_Rgenes.RData"))

library(extrafont)

cognates <- list("AVR_Pia" = c("Pi_a"),
                 "AVR_ii" = c("Pi_i"),
                 "ACE1" = c("pi_d"),
                 "AVR_Pib" = c("Pi_b"),
                 "AVR_Pizt" = c("Pi_z"), # this one
                 "AVR_Pita1" = c("Pi_ta_Pi_ta2_Ptr"),
                 "AVR_Piks" = c("Pi_km", "Pi_ks"))


# fix the year labels
Rgenes_AVR_over_time.lag <- lapply(Rgenes_AVR_over_time.lag,
                                   \(x) paste0("(",
                                               substr(x, 1, 2),
                                               "-",
                                               substr(x, 7,8),
                                               " y.)"))


#### Plot Figure 2 first, which is DROP Pi9, and exclude Pi_km and Pi_ks

set.seed(123)

f2.A <- ggplot2::ggplot(AVR.Global |>
                  dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
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
                      ncol = 1) +
  ggplot2::geom_line(data = AVR.split.Global |>
                       dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.Global |>
                         dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
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
                 panel.spacing.x = ggplot2::unit(1, "lines"),
                 strip.text = ggplot2::element_text(face = "italic"),
                 plot.title = ggplot2::element_text(hjust = 0.5),
                 axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 12, b = 0, l = 0))) +
  ggplot2::scale_y_continuous(name = "Predicted Prevalence of AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25)) +
  ggplot2::scale_x_continuous(name = "Year World Isolate\nSample Collected",
                              limits = c(1970, 2020)) +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             # values = viridisLite::plasma(2),
                             values = c("purple", "lightgrey"),
                             labels = c("Absent", "Present"))


f2.B <- ggplot2::ggplot(AVR.USA |>
                          dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
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
                      ncol = 1) +
  ggplot2::geom_line(data = AVR.split.USA |>
                       dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.USA |>
                         dplyr::filter(AVR %in% c("AVR-Pi9", "AVR-Pib")),
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
                 panel.spacing.x = ggplot2::unit(1, "lines"),
                 strip.text = ggplot2::element_text(face = "italic"),
                 plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::scale_y_continuous(name = "Predicted Prevalence of AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25)) +
  ggplot2::xlab("Year U.S. Isolate\nSample Collected") +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             # values = viridisLite::plasma(2),
                             values = c("purple", "lightgrey"),
                             labels = c("Absent", "Present"))

f2.C <- ggplot2::ggplot(data = Rgenes |>
                          dplyr::filter(R %in% c("Pi-b")),
                       ggplot2::aes(x = Year,
                                    y = proportion,
                                    group = State)) +
  ggplot2::geom_point(ggplot2::aes(fill = State),
                      color = "black", alpha = 0.5, pch = 21, size = 0.75) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 legend.position.inside = c(0.84, 0.85),
                 panel.spacing.x = ggplot2::unit(2, "lines"),
                 strip.text = ggplot2::element_text(face = "italic"),
                 axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 12))) +
  ggplot2::geom_line(data = Rgenes.loess |>
                       dplyr::filter(R %in% c("Pi-b")),
                     linewidth = 0.2,
                     ggplot2::aes(x = Year,
                                  y = pred,
                                  group = State,
                                  color = State)) +
  ggplot2::scale_y_continuous(name = "Proportion of Rice Acres with R Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25),
                              position = "right") +
  ggplot2::scale_x_continuous(name = "Year Rice Planted\nBy U.S. State", limits = c(1960L, 2020L)) +
  ggplot2::facet_wrap(~R, ncol = 1) +
  ggplot2::geom_ribbon(data = Rgenes.loess |>
                         dplyr::filter(R %in% c("Pi-b")),
                       ggplot2::aes(x = Year,
                                    ymin = low,
                                    ymax = high,
                                    group = State,
                                    color = State,
                                    fill = State),
                       alpha = 0.2,
                       inherit.aes = FALSE,
                       # color = "black",
                       linewidth = 0.2)

f2.AB <- ggpubr::ggarrange(f2.A, f2.B + ggpubr::rremove("ylab") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                           common.legend = TRUE, legend = "none",
                           widths = c(1.1, 1),
                           align = "h",
                           labels = c("A", "B"),
                           font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))


f2.ABC <- ggpubr::ggarrange(f2.AB, f2.C,
                            common.legend = FALSE,
                            # align = "h",
                            widths = c(2, 1),
                            labels = c("", "C"),
                            font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

ggplot2::ggsave(plot = f2.ABC, filename = here::here("figures", "Figure_S1.png"),
                dpi = 600, width = 5, height = 3, scale = 1.75)


