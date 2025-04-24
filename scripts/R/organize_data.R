here::i_am("scripts/R/organize_data.R")

dataset.S1 <- readxl::read_xlsx(here::here("inputs", "NP Dataset S1.xlsx"), range = "A6:M846",
                                na = c("NA", "-", ""))
dataset.S2.A <- readxl::read_xlsx(here::here("inputs", "NP Dataset S2.xlsx"), range = "A7:D1778",
                                  na = c("NA", "NAc"))
dataset.S2.B <- readxl::read_xlsx(here::here("inputs", "NP Dataset S2.xlsx"), range = "J7:L179")

dataset.S3 <- readxl::read_xlsx(here::here("inputs", "NP Dataset S3.xlsx"), range = "A7:S400",
                                .name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE))
table.S2 <- readxl::read_xlsx(here::here("inputs", "Table S2.xlsx"))

# set column name references to NULL to avoid warnings
Varieties <- NULL
Year <- NULL
acres <- NULL
State <- NULL
HarvestedAcres <- NULL
Variety <- NULL
pred <- NULL
`Isolate name` <- NULL

dataset.S2.A <- dataset.S2.A |>
  dplyr::rename("HarvestedAcres" = "Harvested Acresa") |>
  dplyr::filter(!is.na(HarvestedAcres))

# fix some mismatched case names
dataset.S2.A$Varieties[dataset.S2.A$Varieties == "Bluebonnet"] <- "BlueBonnet"
dataset.S2.A$Varieties[dataset.S2.A$Varieties == "LaBelle"] <- "Labelle"

# get all years, state, and variety combinations
combns <- expand.grid("Varieties" = unique(dataset.S2.A$Varieties),
                      "HarvestedAcres" = NA_integer_,
                      "Year" = 1959:2020,
                      "State" = c("AR", "LA", "MS", "TX"))

# attach the missing values
combns_keep <- combns[!apply(combns[, -c(2)], 1, paste0, collapse = "") %in%
                           apply(dataset.S2.A[, -c(2)], 1, paste0, collapse = ""), ]
dataset.S2.A <- rbind(dataset.S2.A, combns_keep) |>
  dplyr::arrange(Year, State, Varieties)

rm(combns, combns_keep)

dataset.S2.A.split <- split(dataset.S2.A[, c(1:3)], dataset.S2.A$State) |>
  lapply(\(x) {x2 <- x |>
           tidyr::pivot_wider(names_from = Year, values_from = HarvestedAcres);
         nm <- x2$Varieties;
         x2 <- as.matrix(x2[, -c(1)]);
         # set varieties that are all zero to just zero
         x2[rowSums(x2, na.rm = TRUE) == 0, ] <- 0L;
         x2[, 1] <- 0L;
         x2[, ncol(x2)] <- 0L;
         # set 3 years prior to the first occurrence to zero
         left <- data.frame("row" = 1:nrow(x2),
                            "col" = pmax(apply(x2 != 0, 1, which.max) - 3L, 1)) |>
           as.matrix();
         right <- data.frame("row" = 1:nrow(x2),
                             "col" = ncol(x2) - pmax(apply(x2[, ncol(x2):1] != 0, 1, which.max) - 3L, 1) + 1L) |>
           as.matrix();
         x2[left] <- 0L;
         x2[right] <- 0L;
         rownames(x2) <- nm;
         return(x2)})

# for LA, MS, and TX, replace 1977-1981 with 1976
# and replace 1985-1989 with 1990
dataset.S2.A.split$LA[, as.character(1977:1981)] <- dataset.S2.A.split$LA[, as.character(1976)]
dataset.S2.A.split$LA[, as.character(1985:1989)] <- dataset.S2.A.split$LA[, as.character(1990)]

dataset.S2.A.split$MS[, as.character(1977:1981)] <- dataset.S2.A.split$MS[, as.character(1976)]
dataset.S2.A.split$MS[, as.character(1985:1989)] <- dataset.S2.A.split$MS[, as.character(1990)]

dataset.S2.A.split$TX[, as.character(1977:1981)] <- dataset.S2.A.split$TX[, as.character(1976)]
dataset.S2.A.split$TX[, as.character(1985:1989)] <- dataset.S2.A.split$TX[, as.character(1990)]

# do some linear interpolation on the acreages
xat <- colnames(dataset.S2.A.split$AR)

dataset.S2.A.split <- lapply(dataset.S2.A.split,
                    \(z) apply(z, 1, \(r) approx(x = xat,
                                      y = r,
                                      xout = 1960:2018, rule = 2)$y) |>
           t()) |>
  lapply(\(z) {colnames(z) <- 1960:2018 ; return(z)}) |>
  lapply(floor)

dataset.S2.A <- lapply(dataset.S2.A.split,
       \(x) as.data.frame(x) |>
         tibble::rownames_to_column("Varieties") |>
         tidyr::pivot_longer(-Varieties, names_to = "Year", values_to = "HarvestedAcres")) |>
  purrr::list_rbind(names_to = "State") |>
  dplyr::select(Varieties, HarvestedAcres, Year, State)

rm(dataset.S2.A.split)

# and remove some unnecessary entries that have R genes recorded
# but do not have acreage
table.S2 <- subset(table.S2, !Cultivars %in% setdiff(table.S2$Cultivars, dataset.S2.A$Varieties))

# conver to 0/1 presence absence matrix
rice <- list()
rice$Rgenes <- as.matrix(table.S2[, 2:9] == "1")
rice$Rgenes[] <- as.integer(rice$Rgenes)
rownames(rice$Rgenes) <- table.S2$Cultivars

# set all blanks to zero, but only if the row is not "unpublished"
rice$Rgenes[which(is.na(rice$Rgenes), arr.ind = TRUE)] <- 0L
rice$Rgenes[grepl("unpublished", table.S2$`References and notes`), ] <- NA_integer_
rice$Rgenes <- rice$Rgenes[sort(rownames(rice$Rgenes)), ]

# calculate acreages
# make a wide df where years are cols and varieties are rows
temp.acreage <- dataset.S2.A |>
  dplyr::group_by(Varieties, Year) |>
  dplyr::summarize(acres = sum(HarvestedAcres, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(Year) |>
  tidyr::pivot_wider(names_from = Year, values_from = acres, values_fill = 0) |>
  dplyr::arrange(Varieties)

# convert to a matrix and store in the data list
temp.acreage.mat <- data.matrix(temp.acreage[, -c(1)])
rownames(temp.acreage.mat) <- temp.acreage$Varieties
rice$acreage <- temp.acreage.mat
rm(temp.acreage, temp.acreage.mat)

#####

# calculate the acreage in cultivars with, without, and unknown
# for the R genes
acreage.sus <- t(ifelse(!is.na(rice$Rgenes) & rice$Rgenes == 0, 1, 0)) %*% rice$acreage
acreage.res <- t(ifelse(!is.na(rice$Rgenes) & rice$Rgenes == 1, 1, 0)) %*% rice$acreage
acreage.unk <- t(ifelse(is.na(rice$Rgenes), 1, 0)) %*% rice$acreage
acreage.tot <- acreage.sus + acreage.res + acreage.unk
rice$Racreage <- list("susceptible" = acreage.sus,
                      "resistant" = acreage.res,
                      "unknown" = acreage.unk)

# calculate proportions for each, keeping the unknowns in the counts
rice$Racreage.prop <- lapply(rice$Racreage,
                             \(x) x/acreage.tot)

# calculate proportions for +/- R gene, excluding unknowns
rice$Racreage.prop.known <- lapply(rice$Racreage[1:2],
                                   \(x) x/(acreage.sus + acreage.res))

rm(acreage.sus, acreage.res, acreage.unk, acreage.tot)

# calculate additional quantities
temp.acreage.list <- dataset.S2.A |>
  split(dataset.S2.A$State) |>
  lapply(\(x) {
    temp.acreage <- x |>
      dplyr::group_by(Varieties, Year) |>
      dplyr::summarize(acres = sum(HarvestedAcres, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(Year) |>
      tidyr::pivot_wider(names_from = Year, values_from = acres, values_fill = 0) |>
      dplyr::arrange(Varieties)
    temp.acreage.mat <- data.matrix(temp.acreage[, -c(1)])
    rownames(temp.acreage.mat) <- temp.acreage$Varieties
    temp.acreage.out <- rice$acreage
    temp.acreage.out[] <- 0
    temp.acreage.out[rownames(temp.acreage.mat), colnames(temp.acreage.mat)] <- temp.acreage.mat
    return(temp.acreage.out)
  })

rice$acreage.bystate <- temp.acreage.list

rm(temp.acreage.list)

acreage.sus <- lapply(rice$acreage.bystate,
                      \(x) t(ifelse(!is.na(rice$Rgenes) & rice$Rgenes == 0, 1, 0)) %*% x)
acreage.res <- lapply(rice$acreage.bystate,
                      \(x) t(ifelse(!is.na(rice$Rgenes) & rice$Rgenes == 1, 1, 0)) %*% x)
acreage.unk <- lapply(rice$acreage.bystate,
                      \(x) t(ifelse(is.na(rice$Rgenes), 1, 0)) %*% x)

acreage.tot <- purrr:::pmap(list(acreage.sus, acreage.res, acreage.unk),
                            \(x,y,z) x+y+z)

rice$Racreage.bystate <- purrr:::pmap(list(acreage.sus, acreage.res, acreage.unk),
                                      \(x,y,z) list("susceptible" = x,
                                                    "resistant" = y,
                                                    "unknown" = z))

rm(acreage.sus, acreage.res, acreage.unk, acreage.tot)


rice$Racreage.prop.bystate <- lapply(rice$Racreage.bystate,
                                     \(y) lapply(y, \(x) x/(y$susceptible + y$resistant + y$unknown)))

rice$Racreage.prop.known.bystate <- lapply(rice$Racreage.bystate,
                                           \(y) lapply(y[1:2], \(x) x/(y$susceptible + y$resistant)))


### create the dataset for the pathogen
# set
# FL -> TX
# MO -> AR
# TN -> AR

pathogens.SSR <- dataset.S3 |>
  dplyr::rename("Isolate name" = "...1") |>
  dplyr::select(`Isolate name`, dplyr::starts_with("pyrm")) |>
  dplyr::mutate_all(\(x) ifelse(x == "NA", NA_character_, x)) |>
  dplyr::mutate_at(-1, as.integer) |>
  as.data.frame()

# 6 obs excluded because they are before 1974
# 11 excluded because they are missing state information

pathogens <- dataset.S1 |>
  dplyr::select(1,2,3,4, 6:13) |>
  dplyr::filter(Year >= 1970, Year <= 2018, !is.na(State)) |>
  dplyr::mutate(State = stringr::str_replace_all(State,
                                                 "FL", "TX")) |>
  dplyr::mutate(State = stringr::str_replace_all(State,
                                                 "MO", "AR")) |>
  dplyr::mutate(State = stringr::str_replace_all(State,
                                                 "TN", "AR")) |>
  dplyr::filter(State %in% c("AR", "LA", "MS", "TX")) |>
  as.data.frame() |>
  dplyr::left_join(pathogens.SSR, by = "Isolate name")

rm(pathogens.SSR)
# append the SSRs, if available


xat <- colnames(rice$Racreage.prop.known.bystate$AR$resistant)

rice$Racreage.prop.known.bystate.interpol <- lapply(rice$Racreage.prop.known.bystate,
                                            \(z1) lapply(z1,
                             \(z) apply(z, 1, \(r) approx(x = xat,
                                                          y = r,
                                                          xout = 1960:2018, rule = 2)$y) |>
                               t()) |>
  lapply(\(z) {colnames(z) <- 1960:2018 ; return(z)}))

saveRDS(rice, file = here::here("inputs", "rice.RDS"))
saveRDS(pathogens, file = here::here("inputs", "pathogen.RDS"))

