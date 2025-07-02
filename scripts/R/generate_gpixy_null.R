
# set here:: ------------------------------------------------------------------------------------------------------

here::i_am("scripts/R/generate_gpixy_null.R")

# set ignores -----------------------------------------------------------------------------------------------------

CHR <- NULL
START <- NULL
STOP <- NULL
label <- exon <- NAME <- gene <- NULL

# only look at exons
regions <- read.table(here::here("inputs", "exons.bed.gz"),
             sep = "\t") |>
  setNames(c("CHR", "START", "STOP", "NAME", "SCORE", "STRAND")) |>
  tidyr::separate(NAME,
                  into = c("label", "gene", "exon"),
                  remove = FALSE) |>
  dplyr::mutate("gene" = sprintf("%05d", as.integer(gene))) |>
  dplyr::select(-label, -exon) |>
  dplyr::mutate("REGION" = paste0(CHR, ":", START + 1L, "-", STOP)) |>
  dplyr::mutate("pi" = NA_real_,
                "khat" = NA_real_,
                "S" = NA_real_,
                "N" = NA_real_,
                "M"  = NA_real_,
                "d" = NA_real_,
                "D" = NA_real_)

regions.results.template <- regions |>
  dplyr::select(gene, 9:15) |>
  dplyr::distinct()

nrep.gene <- 100L

regions.split   <- split(regions[, c("CHR", "START", "STOP")], regions$gene)
regions.results <- split(regions.results.template[, -1], regions.results.template$gene) |>
  lapply(\(x) x[rep(1, nrep.gene), ])

rm(regions.results.template)

# load in the data ------------------------------------------------------------------------------------------------

PATHTODATA <- "../"

data <- data.table::fread(fs::path(PATHTODATA, "All.txt.gz")) |>
  as.matrix()

map <- data.table::fread(fs::path(PATHTODATA, "All.map.txt.gz"))


# load in the sample metadata -------------------------------------------------------------------------------------
meta <- read.table(here::here("inputs", "meta_combined.txt"), header = TRUE,
                   comment.char = "", sep = "\t") |>
  dplyr::mutate(Domestic = ifelse(Country == "USA", TRUE, FALSE))

if (gene == "Pita1") {
  # size 22
  meta$Group <- ifelse(meta$Year < 1995,
                       "Pre 1995",
                       ifelse(meta$Year <= 2005,
                              "1995 - 2005",
                              "Post 2005")) |>
    factor(levels = c("Pre 1995", "1995 - 2005", "Post 2005"))
} else if (gene == "Piks") {
  # size 27
  meta$Group <- ifelse(meta$Year < 1990,
                       "Pre 1990",
                       ifelse(meta$Year <= 2010,
                              "1990 - 2010",
                              "Post 2010")) |>
    factor(levels = c("Pre 1990", "1990 - 2010", "Post 2010"))
}