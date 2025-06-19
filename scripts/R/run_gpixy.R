here::i_am("scripts/R/run_gpixy.R")

CHR <- NULL
START <- NULL
STOP <- NULL

regions <- rbind(
  read.table(here::here("inputs", "exons.bed.gz"),
             sep = "\t"),
  read.table(here::here("inputs", "introns.bed.gz"),
             sep = "\t"),
  read.table(here::here("inputs", "intergenic.bed.gz"),
             sep = "\t")) |>
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

PATHTODATA <- "../../"

data <- data.table::fread(fs::path(PATHTODATA, "All.txt.gz")) |>
  as.matrix()

map <- data.table::fread(fs::path(PATHTODATA, "All.map.txt.gz"))

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

# setwd("../scripts/R/genomicpp")
# devtools::document()
# devtools::load_all()

maxit = length(strides$NAME)

regions <- regions[match(strides$NAME, regions$NAME), ]

pb = txtProgressBar(min = 0, max = maxit, initial = 0,
                    width = 80, style = 3)
i = 0

system.time(
for (count in c(1:length(strides$NAME))[1:maxit]) {
  idx1 <- strides[count, "START"]
  idx2 <- strides[count, "STOP"]

  x.mat <- collapse::ss(data, i = idx1:idx2)

  regions[count, 8:14] <- genomicpp::rcpp_parallel_tajimas_d(x.mat)

  i = i + 1

  setTxtProgressBar(pb,i)
}
)

close(pb)

options(scipen = 999)

write.table(regions,
            file = here::here("outputs", "interval_popgen.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

# now run on the aligned sequences...
