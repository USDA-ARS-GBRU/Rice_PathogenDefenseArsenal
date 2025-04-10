regions <- rbind(
  read.table("/project/90daydata/gbru_fy17_005/Yulin/intervals/exons.bed",
             sep = "\t"),
  read.table("/project/90daydata/gbru_fy17_005/Yulin/intervals/introns.bed",
             sep = "\t"),
  read.table("/project/90daydata/gbru_fy17_005/Yulin/intervals/intergenic.bed",
             sep = "\t"))|>
  setNames(c("CHR", "START", "STOP", "NAME", "SCORE", "STRAND")) |>
  dplyr::mutate("REGION" = paste0(CHR, ":", START + 1L, "-", STOP)) |>
  dplyr::mutate("pi" = NA_real_,
                "khat" = NA_real_,
                "S" = NA_real_,
                "N" = NA_real_,
                "M"  = NA_real_,
                "d" = NA_real_,
                "D" = NA_real_)

regions <- regions[gtools::mixedorder(regions$REGION), ]

data <- data.table::fread("/project/90daydata/gbru_fy17_005/Yulin/pixy_new/All.txt.gz") |>
  as.matrix()

map <- data.table::fread("/project/90daydata/gbru_fy17_005/Yulin/pixy_new/All.map.txt.gz")

colnames(map) <- c("CHR", "POS")
map$POS <- map$POS - 1L
map$idx <- 1:nrow(map)

map.split <- split(map[, -c(1)], map$CHR)
regions.split       <- split(regions[, c("START", "STOP", "NAME")], regions$CHR)

strides <- purrr::map2(map.split,
                       regions.split,
                       \(m, r) data.frame("START" = m$idx[match(r$START, m$POS)],
                                          "STOP"  = m$idx[match(r$STOP - 1L, m$POS)],
                                          "NAME" = r$NAME)) |>
  lapply(\(x) x[complete.cases(x), ]) |>
  purrr::list_rbind()

# require a segment to be at least 5 bases
strides <- subset(strides, (STOP - START + 1L) >= 5)

setwd("~/rcpp_examples")
# devtools::document()
# devtools::load_all()

maxit = length(strides$NAME)

regions <- regions[match(strides$NAME, regions$NAME), ]

pb = txtProgressBar(min = 0, max = maxit, initial = 0,
                    width = 80, style = 3)
i = 0


# I think we can make this so much faster if we don't use bcftools...
system.time(
for (count in c(1:length(strides$NAME))[1:maxit]) {
  idx1 <- strides[count, "START"]
  idx2 <- strides[count, "STOP"]

  # x.mat <- data[idx1:idx2, ]
  x.mat <- collapse::ss(data, i = idx1:idx2)

  regions[count, 8:14] <- genomicpp::rcpp_parallel_tajimas_d(x.mat)

  i = i + 1

  setTxtProgressBar(pb,i)
}
)

close(pb)

options(scipen=999)

write.table(regions,
            file = "/project/90daydata/gbru_fy17_005/Yulin/grantpixy.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

regions <- read.table("/project/90daydata/gbru_fy17_005/Yulin/grantpixy.txt",
                      header = TRUE, sep = "\t")
# look at our genes:
regions |>
  dplyr::filter(grepl("_4881_", NAME) |
                  grepl("_5074_", NAME) |
                  grepl("_10623_", NAME) |
                  grepl("_11753_", NAME) |
                  grepl("_12522_", NAME) |
                  grepl("_12603_", NAME)) |> View()

quantile(regions$D, 0:20/20)


ecdf(log10(ifelse(regions$pi == 0,
                  min(regions$pi[regions$pi > 0])/10,
                  regions$pi))) |>
  plot(xlim = c(-8, 0), verticals = TRUE)


ecdf(log10(ifelse(regions$pi == 0,
                  min(regions$pi[regions$pi > 0])/10,
                  regions$pi))) |>
  plot(xlim = c(-8, 0),
       verticals = TRUE)

regions <- regions |>
  tidyr::separate(NAME, into = c("type", "extra"), extra = "merge",
                  remove = FALSE)

ggplot2::ggplot(data = subset(regions, D < 5),
                ggplot2::aes(x = START/1e6,
                             y = D,
                             fill = type)) +
  ggplot2::geom_point(alpha = 0.5, color = "black", shape = 21) +
  ggplot2::geom_smooth(ggplot2::aes(color = type,
                                    linetype = type)) +
  ggplot2::facet_wrap(~CHR, nrow = 2, scale = "free_x") +
  ggplot2::theme_minimal() +
  ggplot2::geom_abline(intercept = mean(regions$D),
                       slope = 0,
                       col = "red", linetype = "dashed")



ggplot2::ggplot(data = subset(regions, D < 5),
                ggplot2::aes(x = START/1e6,
                             y = pi,
                             fill = type)) +
  ggplot2::geom_point(alpha = 0.5, color = "black", shape = 21) +
  ggplot2::geom_smooth(ggplot2::aes(color = type,
                                    linetype = type)) +
  ggplot2::facet_wrap(~CHR, nrow = 2, scale = "free_x") +
  ggplot2::theme_minimal() +
  ggplot2::geom_abline(intercept = mean(regions$D),
                       slope = 0,
                       col = "red", linetype = "dashed") +
  ggplot2::scale_y_continuous(trans=scales::pseudo_log_trans(base = 10,
                                                             sigma = 0.0001),
                              limits = c(0, 0.15),
                              breaks = c(0, 0.0001, 0.0005, 0.001,
                                         0.005, 0.01, 0.02, 0.03, 0.04,
                                         0.06, 0.08, 0.12),
                              labels = scales::label_number(drop0trailing=TRUE))
