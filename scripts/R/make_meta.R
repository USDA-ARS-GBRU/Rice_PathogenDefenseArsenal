here::i_am("scripts/R/make_meta.R")

# samples <- read.csv("Documents/projects/Rice_PathogenDefenseArsenal_Working/information_needed2.csv")
# x <- lapply(samples$GenBank, gbfetch::fetch_metadata)
# xl <- do.call(rbind, x)

meta <- read.table(here::here("sequences", "meta.txt"), header = FALSE,
                   comment.char = "", col.names = "Sample")

Pita1 <- read.table(here::here("sequences", "AVR_Trees", "Pita1", "seqnames.txt"),
                   header = FALSE, comment.char = "")$V1
Piks <- read.table(here::here("sequences", "AVR_Trees", "Piks", "seqnames.txt"),
                   header = FALSE, comment.char = "")$V1
meta$Pita <- meta$Sample %in% Pita1
meta$Piks <- meta$Sample %in% Piks

SRR.meta <- read.table(here::here("inputs", "SRR_meta.txt"), header = TRUE, sep = "\t",
                       check.names = FALSE, col.names = c("PRJNA", "Sample", "Year", "Country"))[, -c(1)] |>
  tibble::add_column("State" = NA_character_, .after = 3)


Wang.meta1 <- readRDS(here::here("inputs", "pathogen.RDS"))[, c(1, 3:4)] |>
  setNames(c("Sample", "Year", "State")) |>
  tibble::add_column("Country" = "USA", .before = 3)

Wang.meta2 <- readRDS(here::here("inputs", "pathogen.RDS"))[, 2:4] |>
  setNames(c("Sample", "Year", "State")) |>
  tibble::add_column("Country" = "USA", .before = 3)

Wang.meta <- rbind(Wang.meta1, Wang.meta2) |>
  dplyr::distinct()

meta2 <- rbind(SRR.meta, Wang.meta)

YearXX <- NULL
Year <- NULL
Sample <- NULL
Country <- NULL

meta.more <- read.csv(here::here("..", "Rice_PathogenDefenseArsenal_Working",
                                 "Original", "phylo", "rericeblastpaperdiscussion",
                                 "US AVRPita isolates with collection period.csv")) |>
  setNames(c("Sample", "Country", "YearXX")) |>
  dplyr::mutate(Year = as.integer(substr(YearXX, 1, 4))) |>
  dplyr::mutate(Year = ifelse(substr(YearXX, 6,6) == "e",
                              Year + 2L,
                              Year + 7L)) |>
  dplyr::select(Sample, Year, Country) |>
  dplyr::mutate(State = NA_character_) |>
  dplyr::mutate(Country = ifelse(Country == "US", "USA", Country)) |>
  dplyr::filter(!Sample %in% meta2$Sample)


meta.more2 <- read.csv(here::here("..", "Rice_PathogenDefenseArsenal_Working",
                                 "Original", "phylo", "rericeblastpaperdiscussion",
                                 "US AVRPik collection perios.csv"),
                       header = FALSE) |>
  setNames(c("YearXX", "Sample")) |>
  dplyr::mutate(Sample = trimws(Sample)) |>
  dplyr::mutate(Year = as.integer(substr(YearXX, 1, 4))) |>
  dplyr::mutate(Year = ifelse(substr(YearXX, 6,6) == "e",
                              Year + 2L,
                              Year + 7L)) |>
  dplyr::mutate(Country = "USA") |>
  dplyr::select(Sample, Year, Country) |>
  dplyr::mutate(State = NA_character_) |>
  dplyr::filter(!Sample %in% meta.more$Sample,
                !Sample %in% meta2$Sample)

meta2 <- rbind(meta2, rbind(meta.more, meta.more2))
meta2$Year <- as.integer(meta2$Year)

meta <- dplyr::left_join(meta, meta2)

# update some manually
meta[meta$Sample == "CG4_2008_", "Year"] <- 2008L
meta[meta$Sample == "CG4_2008_", "Country"] <- "USA"

meta[meta$Sample == "14L21_1_REP", "Year"] <- 2014L
meta[meta$Sample == "14L21_1_REP", "Country"] <- "USA"
meta[meta$Sample == "14L21_1_REP", "State"] <- "LA"

dplyr::filter(meta, is.na(Country)) |> clipr::write_clip()

write.table(meta, file = here::here("inputs", "meta_combined.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t",
            quote = FALSE)
