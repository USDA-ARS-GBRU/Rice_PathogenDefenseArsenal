setwd("../../inputs")

rice <- readRDS(file = "rice.RDS")

Year <- NULL
State_R <- NULL
Pred <- NULL

d <- lapply(rice$Racreage.prop.known.bystate.interpol, \(x) x$resistant) |>
  lapply(\(x) data.frame(x, check.names = FALSE) |>
           tibble::rownames_to_column("R")) |>
  purrr::list_rbind(names_to = "State") |>
  tidyr::pivot_longer(-c(1:2), names_to = "Year", values_to = "proportion") |>
  dplyr::mutate(Year = as.integer(Year))

d.loess <- split(d[, c("R", "Year", "proportion")], d$State) |>
  lapply(\(x) split(x[, c("Year", "proportion")], x$R)) |>
  lapply(\(State) lapply(State, \(R) stats::loess(proportion ~ Year, data = R, span = 0.2))) |>
  purrr::list_flatten() |>
  lapply(\(mod) cbind(data.frame("Year" = 1960:2018),
                      data.frame("Pred" = predict(mod, data.frame("Year" = 1960:2018))))) |>
  purrr::list_rbind(names_to = "State_R") |>
  tidyr::separate(State_R, into = c("State", "R"), sep = "_") |>
  dplyr::mutate(Pred = ifelse(Pred < 0, 0, Pred)) |>
  dplyr::mutate(Pred = ifelse(Pred > 1, 1, Pred)) |>
  tibble::remove_rownames()
  
saveRDS(d.loess, file = "rice_smooth.RDS") 

library(extrafont)
extrafont::font_import()

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
                          y = Pred,
                          group = State,
                          color = State)) +
  ggplot2::scale_y_continuous(name = "Proportion of Rice Acres with R Gene",
                              limits = c(0, 1)) +
  ggplot2::scale_x_continuous(limits = c(1960L, 2020L)) +
  ggplot2::facet_wrap(~R, nrow = 2) +
  ggplot2::ggtitle("Rice R Gene Frequency by State")

ggplot2::ggsave(g,
                filename = "../figures/Rgenes_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")
