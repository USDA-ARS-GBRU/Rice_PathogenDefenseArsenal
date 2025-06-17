setwd("../../inputs")

rice <- readRDS(file = "rice.RDS")

Year <- NULL
State_R <- NULL
pred <- NULL
se <- NULL
high <- low <- NULL

d <- lapply(rice$Racreage.prop.known.bystate.interpol, \(x) x$resistant) |>
  lapply(\(x) data.frame(x, check.names = FALSE) |>
           tibble::rownames_to_column("R")) |>
  purrr::list_rbind(names_to = "State") |>
  tidyr::pivot_longer(-c(1:2), names_to = "Year", values_to = "proportion") |>
  dplyr::mutate(Year = as.integer(Year))

d.loess <- split(d[, c("R", "Year", "proportion")], d$State) |>
  lapply(\(x) split(x[, c("Year", "proportion")], x$R)) |>
  lapply(\(State) lapply(State, \(R) stats::loess(car::logit(proportion, adjust = 0.0001) ~ Year, data = R, span = 0.2))) |>
  purrr::list_flatten() |>
  lapply(\(mod) cbind(data.frame("Year" = 1960:2018),
                      data.frame("Pred" = predict(mod, data.frame("Year" = 1960:2018),
                                                  se = TRUE)[1:2]) |>
                        setNames(c("pred", "se")))) |>
  purrr::list_rbind(names_to = "State_R") |>
  tidyr::separate(State_R, into = c("State", "R"), sep = "_") |>
  dplyr::mutate(low = plogis(pred - se),
                high = plogis(pred + se),
                pred = plogis(pred)) |>
  dplyr::mutate(pred = ifelse(pred <= 1e-4,
                              0,
                              ifelse(pred > 1 - 1e-4,
                                     1.0,
                                     pred)),
                high = ifelse(high <= 1e-4,
                              0,
                              ifelse(high > 1 - 1e-4,
                                     1.0,
                                     high)),
                low = ifelse(low <= 1e-4,
                              0,
                              ifelse(low > 1 - 1e-4,
                                     1.0,
                                     low))) |>
  tibble::remove_rownames()

saveRDS(d.loess, file = "rice_smooth.RDS")

library(extrafont)
# extrafont::font_import()

g <- ggplot2::ggplot(data = d,
                ggplot2::aes(x = Year,
                             y = proportion,
                             group = State)) +
  ggplot2::geom_point(ggplot2::aes(fill = State),
                      color = "black", alpha = 0.5, pch = 21, size = 0.75) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 legend.position.inside = c(0.95, 0.85),
                 panel.spacing.x = ggplot2::unit(2, "lines")) +
  ggplot2::geom_line(data = d.loess,
                     linewidth = 1,
             ggplot2::aes(x = Year,
                          y = pred,
                          group = State,
                          color = State)) +
  ggplot2::scale_y_continuous(name = "Proportion of Rice Acres with R Gene",
                              limits = c(0, 1)) +
  ggplot2::scale_x_continuous(limits = c(1960L, 2020L)) +
  ggplot2::facet_wrap(~R, nrow = 2) +
  ggplot2::ggtitle("Rice R Gene Frequency by State") +
  ggplot2::geom_ribbon(data = d.loess,
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

g
ggplot2::ggsave(g,
                filename = "../figures/Rgenes_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")

Rgenes <- d
Rgenes.loess <- d.loess
save(Rgenes, Rgenes.loess, file = "../outputs/source_data_Rgenes.RData")
