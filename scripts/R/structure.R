here::i_am("scripts/R/structure.R")

library(extrafont)
# extrafont::font_import(prompt = FALSE)
extrafont::loadfonts(device="postscript")


ID <- NULL
Decade <- NULL
P <- NULL
AVR <- NULL
probability <- NULL
oldPopulation <- NULL
s1 <- s2 <- s3 <- s4 <- s5 <- NULL
Population <- NULL
Sample <- NULL

pathogen <- readRDS(here::here("inputs", "pathogen.RDS"))

pathogen.SSR <- data.matrix(pathogen[, c(13:22)])
rownames(pathogen.SSR) <- pathogen$`Isolate name`

pathogen.SSR <- pathogen.SSR[rowSums(is.na(pathogen.SSR)) < 10, ]
pathogen.SSR[is.na(pathogen.SSR)] <- -9

pathogen.AVR <- data.matrix(pathogen[, c(5:12)])
rownames(pathogen.AVR) <- pathogen$`Isolate name`
pathogen.AVR <- pathogen.AVR[rownames(pathogen.SSR), ]

write.table(pathogen.SSR,
            file = here::here("inputs", "structure.str"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

pathogen.SSR.meta <- pathogen[match(rownames(pathogen.SSR), pathogen$`Isolate name`), 1:12]


genind <- adegenet::df2genind(pathogen.SSR, ploidy = 1, NA.char = "-9",
                              pop =  pathogen.SSR.meta$State)

set.seed(123)

grp <- adegenet::find.clusters(genind,
                               max.n.clust = 20,
                               n.pca = 1e5,
                               n.iter = 10000,
                               n.start = 100,
                               choose.n.clus = FALSE,
                               criterion = "diffNgroup",
                               pca.center = TRUE,
                               pca.scale = FALSE,
                               n.pca.cores = 1)

grp

dapc1 <- adegenet::dapc(genind,
                        grp$grp,
                        n.pca = 1e5,
                        n.da = 2,
                        var.contrir = TRUE,
                        var.loadings = TRUE,
                        pca.info = TRUE,
                        center = TRUE,
                        scale = FALSE)

adegenet::scatter.dapc(dapc1,
                       col = viridisLite::magma(5),
                       scree.da = TRUE,
                       scree.pca = TRUE,
                       label = NULL,
                       mstree = TRUE,
                       solid = 0.5,
                       pch = 19,
                       cex = 3)

#### copy code from ggcompoplot and modify
da.object <- dapc1
gid <- genind
pal = viridisLite::plasma
cols = 1
posterior <- da.object$posterior

names(dimnames(posterior)) <- c("sample", "population")

to_merge <- data.frame(list(sample = dimnames(posterior)$sample,
                            oldPopulation = adegenet::pop(gid)))

post <- reshape2::melt(posterior, value.name = "probability")
post <- merge(post, to_merge)
if (is.numeric(post$sample)) {
  post$sample <- factor(post$sample, levels = unique(post$sample))
}
if (is.numeric(post$population)) {
  post$population <- factor(post$population, levels = unique(sort(post$population)))
}

# from `ggcompoplot::`
char2pal <- function (x, pal = rainbow)
{
  PAL <- match.fun(pal)
  outPal <- PAL(length(unique(x)))
  names(outPal) <- unique(x)
  return(outPal)
}

if (length(pal) == 1) {
  PAL <- match.fun(pal)
  pal <- char2pal(post$population, PAL)
}

post$Decade <- paste0(floor(pathogen.SSR.meta$Year[match(post$sample, pathogen.SSR.meta$`Isolate name`)] / 10) * 10,
                      "s")
post$Decade <- factor(post$Decade,
                      levels = c("1970s", "1980s", "1990s", "2000s", "2010s"))
# post$Decade <- factor(substr(post$Decade, 3,4),
#                       levels = c("70", "80", "90", "00", "10"))

grp2 <- max.col(dapc1$posterior)
names(grp2) <- rownames(dapc1$posterior)

xord <- dapc1$grp.coord |> dist() |> hclust()
xord <- xord$order
xord <- c(5,1,4,2,3)
grp2[] <- xord[grp2]

post$grp2 <-  grp2[as.character(post$sample)]

post$population <- factor(xord[post$population])

post$sample <- factor(as.character(post$sample),
                      levels = unique(dplyr::arrange(post,
                                                     grp2,
                                                     dplyr::desc(probability))$sample))

g1 <- ggplot2::ggplot(post,
                      ggplot2::aes(x = sample,
                                   fill = population,
                                   y = probability)) +
  ggplot2::geom_bar(stat = "identity",
                    position = "fill",
                    width = 1) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     hjust = 1,
                                                     vjust = 0.5),
                 text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "top",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.3),
                                                           color = "black"),
                 panel.spacing.x = ggplot2::unit(c(0.05,0.05,0.05,0.4,
                                                   0.05,0.4,
                                                   0.05,0.05,0.05,0.4,
                                                   0.05,0.05,0.05,0.05,0.4,
                                                   0.05,0.05,0.05,0.05), "lines"),
                 axis.text.x.bottom = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 strip.text  = ggplot2::element_text(size = 7, family = "Arial", color = "black")) +
  ggplot2::scale_y_continuous(name = "Membership Probability", expand = c(0, 0)) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggh4x::facet_nested(~Decade+grp2,
                      space = "free_x",
                      scales = "free_x",) +
  ggplot2::scale_fill_manual(name = "DAPC Cluster",
                             values = pal)

g1


xord2 <- dist(t(pathogen.AVR)) |> hclust()
xord2 <- rev(xord2$order)
# now make the same plot, but this time with the AVR genes marked
post2 <- post |>
  dplyr::select(sample, Decade, grp2) |>
  dplyr::distinct() |>
  dplyr::left_join(pathogen.AVR |>
                     data.frame(check.names = FALSE) |>
                     tibble::rownames_to_column("sample"),
                   by = "sample") |>
  tidyr::pivot_longer(-c(1:3), names_to = "AVR", values_to = "P") |>
  dplyr::mutate(sample = factor(sample, levels = unique(post$sample))) |>
  dplyr::mutate(P = ifelse(P == 0, "Absent", "Present")) |>
  dplyr::mutate(AVR = factor(AVR,
                             levels = colnames(pathogen.AVR)[xord2]))


g2 <- ggplot2::ggplot(post2 |>
                        dplyr::mutate(AVR = AVR |>
                                        stringr::str_replace_all("Pi_",
                                                                      "Pi") |>
                                        stringr::str_replace_all("_", "-") |>
                                        stringr::str_replace_all("Pita-Pita2-Ptr",
                                                                 "Pi-ta Pi-ta2 Ptr") |>
                                        stringr::str_replace_all("AVR-ii",
                                                                 "AVR-Pii")) |>
                        dplyr::mutate(AVR = factor(AVR,
                                                   levels = c("AVR-Piks",
                                                              "AVR-Pib",
                                                              "AVR-Pi9",
                                                              "AVR-Pizt",
                                                              "AVR-Pita1",
                                                              "AVR-Pii",
                                                              "AVR-Pia",
                                                              "ACE1"))),
                         ggplot2::aes(x = sample,
                                      fill = P,
                                      y = AVR)) +
  ggplot2::geom_raster() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     hjust = 1,
                                                     vjust = 0.5),
                 text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "bottom",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.3),
                                                           color = "black"),
                 panel.spacing.x = ggplot2::unit(c(0.05,0.05,0.05,0.4,
                                                   0.05,0.4,
                                                   0.05,0.05,0.05,0.4,
                                                   0.05,0.05,0.05,0.05,0.4,
                                                   0.05,0.05,0.05,0.05), "lines"),
                 axis.text.x.bottom = ggplot2::element_blank(),
                 strip.background = ggplot2::element_blank(),
                 strip.text = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(face = "italic"),
                 # axis.title.y = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank()) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggh4x::facet_nested(~Decade+grp2,
                      space = "free_x",
                      scales = "free_x") +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             values = c("black", "lightgrey")) +
  ggplot2::ylab("Gene")

g2

gp <- ggpubr::ggarrange(g1, g2, ncol = 1,
                  heights = c(2,1),
                  align = "v",
                  labels = c("A", "B"),
                  font.label = list(size = 16, color = "black", face = "bold", family = "Arial"))

gp
ggplot2::ggsave(gp,
       filename = here::here("figures", "Figure_1.png"),
       width = 7.5,
       height = 3,
       dpi = 600,
       scale = 2,
       units = "in")

post.save <- post |>
  tidyr::pivot_wider(names_from = "population", values_from = "probability", names_prefix = "s") |>
  dplyr::select(sample, Decade, oldPopulation, grp2, s1, s2, s3, s4, s5) |>
  setNames(c("Sample", "Decade", "State", "Most Likely Population",
             "Subpop 1", "Subpop 2", "Subpop 3", "Subpop 4", "Subpop 5")) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("Subpop"), \(x) round(x, 4)))

# merge on the P/A of the genes
post.save.joined <- post.save[, 1:4] |>
  dplyr::rename(Population = 4) |>
  dplyr::mutate(Population = factor(Population)) |>
  dplyr::left_join(post2[, c(1, 2, 4, 5)] |>
                     dplyr::mutate(P = factor(ifelse(P == "Present", 1L, 0L),
                                              levels = c(0L, 1L))) |>
                     tidyr::pivot_wider(names_from = "AVR",
                                        values_from = "P"),
                   by = c("Sample" = "sample",
                          "Decade" = "Decade")) |>
  dplyr::relocate(Population, .before = 2) |>
  dplyr::left_join(pathogen.SSR |>
                     as.data.frame() |>
                     tibble::rownames_to_column("Sample") |>
                     dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) factor(dplyr::na_if(x, -9L)))),
                   by = "Sample")

post.save.joined.supp <- post.save.joined[, c("Sample", "Population", "Decade", "State",
                                              "ACE1", "AVR-Pia", "AVR-ii", "AVR-Pita1",
                                              "AVR-Pizt", "AVR-Pi9", "AVR-Pib", "AVR-Piks",
                                              "pyrm_37",  "pyrm_43", "pyrm_47", "pyrm_63", "pyrm_77",
                                              "pyrm_233", "pyrm_409", "pyrm_427", "pyrm_607", "pyrm_657")] |>
  dplyr::relocate("Population", .after = 4) |>
  dplyr::rename("AVR-Pii" = "AVR-ii",
                "DAPC Population" = "Population") |>
  dplyr::arrange(Sample)

write.table(post.save.joined.supp,
            file = here::here("tables_staged", "Dataset_S1.csv"),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",",
            quote = FALSE,
            na = "-")

# vars <- Filter(\(x) length(unique(x)) > 1L, post.save.joined)
#
# # drop AVR-Pib because the `0` class is way too small
# vars <- subset(vars, select = -`AVR-Pib`)
#
# # rearrange to alphabetize
# vars <- vars[, c(1, 2, 3, 6, 7, 5, 4, 8, 10, 13, 14, 16, 18, 9, 11, 12, 15, 17)]
#
# vars.mat <- matrix(data = NA_real_, nrow = ncol(vars), ncol = ncol(vars),
#                    dimnames = list(colnames(vars), colnames(vars)))
#
# pvals.mat <- matrix(data = NA_real_, nrow = ncol(vars), ncol = ncol(vars),
#                    dimnames = list(colnames(vars), colnames(vars)))
#
# vars.mat[lower.tri(vars.mat)] <- combn(vars,
#                                        2,
#                                        \(x) chisq.test(x = x[[1]],
#                                                        y = x[[2]])$statistic,
#                                        simplify = TRUE)
#
# pvals.mat[lower.tri(pvals.mat)] <- combn(vars,
#                                          2,
#                                          \(x) chisq.test(x = x[[1]],
#                                                          y = x[[2]]) |>
#                                            confintr::cramersv(),
#                                          simplify = TRUE)
#
# H0 <- combn(vars,
#             2,
#             \(x) vapply(1:2000,
#                         \(y) tryCatch(chisq.test(x = sample(x = x[[1]]),
#                                                  y = x[[2]])$statistic,
#                                       error = \(e) NA_real_),
#                         numeric(1)))
#
# xobs <- vars.mat[lower.tri(vars.mat)]
#
# pvals.mat[lower.tri(pvals.mat)] <-  mapply(\(x, y) mean(x < y, na.rm = TRUE),
#                                            as.data.frame(H0),
#                                            xobs) |>
#   unname()
#
# pvals.mat[upper.tri(pvals.mat)] <- t(pvals.mat)[upper.tri(pvals.mat)]
#
# corrplot::corrplot(pvals.mat, is.corr = FALSE, addshade = "positive",
#                    type = "upper", diag = FALSE,
#                    method = "shade",
#                    col = rev(corrplot::COL1("Purples", n = 100)))


# save membership probabilities
post.save |>
  write.csv(file = here::here("outputs", "subpopulations.csv"),
            row.names = FALSE, quote = FALSE)
