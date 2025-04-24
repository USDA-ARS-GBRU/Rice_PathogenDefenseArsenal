here::i_am("scripts/R/structure.R")

library(extrafont)
extrafont::font_import(prompt = FALSE)


ID <- NULL
Decade <- NULL
P <- NULL
AVR <- NULL
probability <- NULL
oldPopulation <- NULL
s1 <- s2 <- s3 <- s4 <- s5 <- NULL

pathogen <- readRDS(here::here("inputs", "pathogen.RDS"))

pathogen.SSR <- data.matrix(pathogen[, c(13:22)])
rownames(pathogen.SSR) <- pathogen$`Isolate name`

pathogen.SSR <- pathogen.SSR[rowSums(is.na(pathogen.SSR)) < 10, ]
pathogen.SSR[is.na(pathogen.SSR)] <- -9

pathogen.AVR <- data.matrix(pathogen[, c(5:12)])
rownames(pathogen.AVR) <- pathogen$`Isolate name`
pathogen.AVR <- pathogen.AVR[rownames(pathogen.SSR), ]


# pathogen.SSR <- pathogen.SSR |> as.data.frame() |>
#   tibble::rownames_to_column("ID") |>
#   dplyr::mutate(ID = stringr::str_replace_all(ID, "/", "")) |>
#   dplyr::mutate(ID = stringr::str_replace_all(ID, " ", ""))
#
# colnames(pathogen.SSR)[1] <- ""

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
                      "'s")
post$Decade <- factor(post$Decade,
                      levels = c("1970's", "1980's", "1990's", "2000's", "2010's"))
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


g2 <- ggplot2::ggplot(post2,
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
                 # axis.title.y = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank()) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggh4x::facet_nested(~Decade+grp2,
                      space = "free_x",
                      scales = "free_x") +
  ggplot2::scale_fill_manual(name = "AVR Status",
                             values = c("black", "lightgrey")) +
  ggplot2::ylab("AVR Gene")

g2

gp <- ggpubr::ggarrange(g1, g2, ncol = 1,
                  heights = c(2,1),
                  align = "v")

ggplot2::ggsave(gp,
       filename = here::here("figures", "AVR_and_structure.png"),
       width = 4,
       height = 2,
       dpi = 600,
       scale = 3,
       units = "in")

# save membership probabilities
post |>
  tidyr::pivot_wider(names_from = "population", values_from = "probability", names_prefix = "s") |>
  dplyr::select(sample, Decade, oldPopulation, grp2, s1, s2, s3, s4, s5) |>
  setNames(c("Sample", "Decade", "State", "Most Likely Population",
             "Subpop 1", "Subpop 2", "Subpop 3", "Subpop 4", "Subpop 5")) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("Subpop"), \(x) round(x, 4))) |>
  write.csv(file = here::here("outputs", "subpopulations.csv"),
            row.names = FALSE, quote = FALSE)
