here::i_am("scripts/R/plot_figure2.R")

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


#### Plot Figure 2 first, which is DROP Pi9, and exclude Pi_km and Pi_ks

set.seed(123)

f2.A <- ggplot2::ggplot(AVR.Global |>
                  dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
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
                       dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.Global|>
                         dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
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
                          dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
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
                       dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
                     ggplot2::aes(x = Year_Collected,
                                  y = pred),
                     color = "red") +
  ggplot2::geom_ribbon(data = AVR.split.USA |>
                         dplyr::filter(! AVR %in% c("AVR-Pi9", "AVR-Piks", "AVR-Pib")),
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
                          dplyr::filter(! R %in% c("Pi-km", "Pi-ks", "Pi-b")) |>
                          dplyr::mutate(R = factor(R, levels = c("pi-d",
                                                                 "Pi-a",
                                                                 "Pi-i",
                                                                 "Pi-ta Pi-ta2 Ptr",
                                                                 "Pi-z"))),
                       ggplot2::aes(x = Year,
                                    y = proportion,
                                    group = State)) +
  ggplot2::geom_point(ggplot2::aes(fill = State),
                      color = "black", alpha = 0.5, pch = 21, size = 0.75) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 legend.position.inside = c(0.86, 0.94),
                 panel.spacing.x = ggplot2::unit(2, "lines"),
                 strip.text = ggplot2::element_text(face = "italic"),
                 axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 12))) +
  ggplot2::geom_line(data = Rgenes.loess |>
                       dplyr::filter(! R %in% c("Pi-km", "Pi-ks", "Pi-b")) |>
                       dplyr::mutate(R = factor(R, levels = c("pi-d",
                                                              "Pi-a",
                                                              "Pi-i",
                                                              "Pi-ta Pi-ta2 Ptr",
                                                              "Pi-z"))),
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
                         dplyr::filter(! R %in% c("Pi-km", "Pi-ks", "Pi-b")) |>
                         dplyr::mutate(R = factor(R, levels = c("pi-d",
                                                                "Pi-a",
                                                                "Pi-i",
                                                                "Pi-ta Pi-ta2 Ptr",
                                                                "Pi-z"))),
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

f2.D <- ggplot2::ggplot(Rgenes_AVR_over_time |>
                          dplyr::filter(AVR != "AVR_Piks",
                                        AVR != "AVR_Pib",
                                        AVR != "AVR_Pi9") |>
                          dplyr::arrange(AVR, R, proportion) |>
                          dplyr::mutate(AVR = factor(AVR,
                                                     levels = c("ACE1",
                                                                "AVR_Pia",
                                                                "AVR_ii",
                                                                "AVR_Pita1",
                                                                "AVR_Pizt"))),
                        ggplot2::aes(x = proportion,
                                     y = prevalence)) +
  ggplot2::geom_line(color = "red") +
  ggplot2::facet_wrap(~AVR,
                      labeller = ggplot2::labeller(AVR = function(x) {
                        x <- as.character(x);
                        return(sapply(x, \(z) paste0(z, " ~ ", cognates[[z]],
                                                     " ", Rgenes_AVR_over_time.lag[[z]])) |>
                                 stringr::str_replace_all("Pi_",
                                                          "Pi") |>
                                 stringr::str_replace_all("_", "-") |>
                                 stringr::str_replace_all("Pita-Pita2-Ptr",
                                                          "Pi-ta Pi-ta2 Ptr") |>
                                 stringr::str_replace_all("AVR-ii",
                                                          "AVR-Pii"));
                      }),
                      ncol = 1) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 panel.spacing.x = ggplot2::unit(2, "lines"),
                 strip.text = ggplot2::element_text(size = 9),
                 axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 12))) +
  ggplot2::scale_y_continuous(name = "Predicted Proportion of Isolates with Cognate AVR Gene",
                              limits = c(-0.1, 1.1),
                              breaks = seq(0, 1, 0.25),
                              position = "right") +
  ggplot2::scale_x_continuous(name = "Prevalence of R Gene\nIn Rice Fields") +
  ggplot2::geom_ribbon(ggplot2::aes(x = proportion,
                                    ymin = pmax(prevalence - SE, 0),
                                    ymax = pmin(prevalence + SE, 1)),
                       alpha = 0.2,
                       inherit.aes = FALSE,
                       color = "black")

f2.AB <- ggpubr::ggarrange(f2.A, f2.B + ggpubr::rremove("ylab") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                           common.legend = TRUE, legend = "none",
                           widths = c(1.1, 1),
                           align = "h",
                           labels = c("A", "B"),
                           font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))


f2.ABC <- ggpubr::ggarrange(f2.AB, f2.C + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                            common.legend = FALSE,
                            align = "h", widths = c(2, 1),
                            labels = c("", "C"),
                            font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

f2.ABCD <- ggpubr::ggarrange(f2.ABC, f2.D, #+ ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"),
                             common.legend = FALSE,
                             align = "h",
                             widths = c(3, 1),
                             labels = c("", "D"),
                             font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

f2.ABCD

ggplot2::ggsave(plot = f2.ABCD, filename = here::here("figures", "Figure_2.png"),
                dpi = 600, width = 7.5, height = 5, scale = 1.75)

State <- Detected <- NULL

## construct Dataset.S2
out <- AVR.Global[, c("ID", "AVR", "Detected", "PRJNA Project", "Year_Collected", "Place/Region")] |>
  tidyr::pivot_wider(names_from = AVR, values_from = Detected) |>
  dplyr::rename("Sample" = "ID")

out2 <- AVR.USA |>
  tidyr::pivot_wider(names_from = AVR, values_from = Detected) |>
  dplyr::mutate("PRJNA Project" = NA_character_, .after = 1) |>
  dplyr::relocate(State, .after = 4) |>
  dplyr::mutate(State = paste0("USA - ", State)) |>
  dplyr::rename("Place/Region" = "State")

out2 <- out2[, colnames(out)]

out <- rbind(out2, out)

rm(out2)

Year_Collected <-  `Place/Region` <- Year <- NULL

write.table(out |>
              dplyr::arrange(Year_Collected, `Place/Region`),
            file = here::here("tables_staged", "Dataset_S2.csv"),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",",
            quote = FALSE,
            na = "-")

rm(out)

out <- Rgenes |>
  tidyr::pivot_wider(names_from = R, values_from = proportion) |>
  dplyr::arrange(State, Year) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.double), \(x) sprintf("%0.4f", x)))

write.table(out,
            file = here::here("tables_staged", "Dataset_S3.csv"),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",",
            quote = FALSE,
            na = "-")

