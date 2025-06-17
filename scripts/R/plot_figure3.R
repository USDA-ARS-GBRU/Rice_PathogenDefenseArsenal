AVR.Global <- AVR.split.Global <- AVR.split.USA <- AVR.USA <- NULL
Rgenes <- Rgenes.loess <-  NULL
Rgenes_AVR_over_time <- Rgenes_AVR_over_time2 <- NULL
Rgenes_AVR_over_time.lag <- Rgenes_AVR_over_time.lag2 <- NULL

AVR <- R <- proportion <- NULL

load(file = "../outputs/source_data_Rgenes_AVR_over_time.RData")
load(file = "../outputs/source_data_AVR_PA_over_time_USA.RData")
load(file = "../outputs/source_data_AVR_PA_over_time_Global.RData")
load(file = "../outputs/source_data_Rgenes.RData")

library(extrafont)
extrafont::loadfonts(device="postscript")

cognates <- list("AVR_Pia" = c("Pi_a"),
                 "AVR_ii" = c("Pi_i"),
                 "ACE1" = c("pi_d"),
                 "AVR_Pib" = c("Pi_b"),
                 "AVR_Pizt" = c("Pi_z"),
                 "AVR_Pita1" = c("Pi_ta_Pi_ta2_Ptr"),
                 "AVR_Piks" = c("Pi_km", "Pi_ks"))


# fix the year labels
Rgenes_AVR_over_time.lag <- lapply(Rgenes_AVR_over_time.lag,
                                   \(x) paste0("(",
                                               substr(x, 1, 2),
                                               "-",
                                               substr(x, 7,8),
                                               " y.)"))

Rgenes_AVR_over_time.lag2 <- sapply(Rgenes_AVR_over_time.lag2,
                                    \(x) paste0("(",
                                                substr(x, 1, 2),
                                                "-",
                                                substr(x, 7,8),
                                                " y.)"))



#### Plot Figure 2 first, which is DROP Pi9, and exclude Pi_km and Pi_ks

set.seed(123)

f3.A <- ggplot2::ggplot(AVR.Global |>
                  dplyr::filter(AVR %in% c("AVR-Piks")),
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
  ggplot2::geom_line(data = AVR.split.Global|>
                       dplyr::filter(AVR %in% c("AVR-Piks")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.Global|>
                         dplyr::filter(AVR %in% c("AVR-Piks")),
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
  ggplot2::scale_y_continuous(name = "Predicted Prevalence\nof AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25)) +
  ggplot2::scale_x_continuous(name = "Year World Isolate\nSample Collected",
                              limits = c(1970, 2020)) +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             # values = viridisLite::plasma(2),
                             values = c("purple", "lightgrey"),
                             labels = c("Absent", "Present"))


f3.B <- ggplot2::ggplot(AVR.USA |>
                          dplyr::filter(AVR %in% c("AVR-Piks")),
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
                       dplyr::filter(AVR %in% c("AVR-Piks")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.USA |>
                         dplyr::filter(AVR %in% c("AVR-Piks")),
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
  ggplot2::scale_y_continuous(name = "Predicted Prevalence\nof AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25)) +
  ggplot2::xlab("Year U.S. Isolate\nSample Collected") +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             # values = viridisLite::plasma(2),
                             values = c("purple", "lightgrey"),
                             labels = c("Absent", "Present"))

f3.C <- ggplot2::ggplot(data = Rgenes |>
                          dplyr::filter(R %in% c("Pi-km", "Pi-ks")),
                       ggplot2::aes(x = Year,
                                    y = proportion,
                                    group = State)) +
  ggplot2::geom_point(ggplot2::aes(fill = State),
                      color = "black", alpha = 0.5, pch = 21, size = 0.75) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 legend.position.inside = c(0.5, 0.88),
                 panel.spacing.x = ggplot2::unit(2, "lines"),
                 strip.text = ggplot2::element_text(face = "italic"),
                 axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 12))) +
  ggplot2::geom_line(data = Rgenes.loess |>
                       dplyr::filter(R %in% c("Pi-km", "Pi-ks")),
                     linewidth = 0.2,
                     ggplot2::aes(x = Year,
                                  y = pred,
                                  group = State,
                                  color = State)) +
  ggplot2::scale_y_continuous(name = "Proportion of Rice Acres\nwith R Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25),
                              position = "left") +
  ggplot2::scale_x_continuous(name = "Year Rice Planted\nBy U.S. State", limits = c(1960L, 2020L)) +
  ggplot2::facet_wrap(~R, nrow = 1) +
  ggplot2::geom_ribbon(data = Rgenes.loess |>
                         dplyr::filter(R %in% c("Pi-km", "Pi-ks")),
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

f3.D <- ggplot2::ggplot(Rgenes_AVR_over_time2,
                          ggplot2::aes(x = Pi_km,
                                       y = AVR_Piks)) +
  ggplot2::geom_line(color = "red") +
  ggplot2::theme_minimal() +
  ggplot2::facet_wrap(~Pi_ks,
                      labeller = ggplot2::labeller(Pi_ks = function(x) {
                        x <- as.character(x);
                        return(paste0("AVR_Piks ~ Pi_km | Pi-ks: ", x, "\n", Rgenes_AVR_over_time.lag2) |>
                                 stringr::str_replace_all("Pi_",
                                                          "Pi") |>
                                 stringr::str_replace_all("_", "-") |>
                                 stringr::str_replace_all("Pita-Pita2-Ptr",
                                                          "Pi-ta Pi-ta2 Ptr") |>
                                 stringr::str_replace_all("AVR-ii",
                                                          "AVR-Pii"));
                      }),
                      ncol = 3) +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 panel.spacing.x = ggplot2::unit(2, "lines"),
                 strip.text = ggplot2::element_text(size = 12),
                 axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 12))) +
  ggplot2::scale_y_continuous(name = "Predicted Proportion of Isolates with Cognate AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25),
                              position = "right") +
  ggplot2::scale_x_continuous(name = "Prevalence of R Gene\nIn Rice Fields") +
  ggplot2::geom_ribbon(ggplot2::aes(x = Pi_km,
                                    ymin = pmax(AVR_Piks - SE, 0),
                                    ymax = pmin(AVR_Piks + SE, 1)),
                       alpha = 0.2,
                       inherit.aes = FALSE,
                       color = "black")

f3.AB <- ggpubr::ggarrange(f3.A, f3.B + ggpubr::rremove("ylab") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                           ncol = 2,
                           common.legend = TRUE, legend = "none",
                           widths = c(1.1, 1),
                           align = "h",
                           labels = c("A", "B"),
                           font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

f3.AB

f3.ABC <- ggpubr::ggarrange(f3.AB, f3.C,
                            common.legend = FALSE,
                            nrow = 2,
                            labels = c("", "C"),
                            font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

f3.ABC

f3.ABCD <- ggpubr::ggarrange(f3.ABC, f3.D, #+ ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                             common.legend = FALSE,
                             align = "h",
                             widths = c(2, 3),
                             labels = c("", "D"),
                             font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

f3.ABCD

ggplot2::ggsave(plot = f3.ABCD, filename = "../figures/Figure_3.png",
                dpi = 600, width = 7.5, height = 3.5, scale = 2)


