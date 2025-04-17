setwd("../../inputs")

ID <- NULL

pathogen <- readRDS("pathogen.RDS")

pathogen.AVR <- data.matrix(pathogen[, c(13:22)])
rownames(pathogen.AVR) <- pathogen$`Isolate name` #|>
# stringr::str_replace_all("/", "")
pathogen.AVR <- pathogen.AVR[rowSums(is.na(pathogen.AVR)) < 10, ]
pathogen.AVR[is.na(pathogen.AVR)] <- -9

# pathogen.AVR <- pathogen.AVR |> as.data.frame() |>
#   tibble::rownames_to_column("ID") |>
#   dplyr::mutate(ID = stringr::str_replace_all(ID, "/", "")) |>
#   dplyr::mutate(ID = stringr::str_replace_all(ID, " ", ""))
# 
# colnames(pathogen.AVR)[1] <- ""

write.table(pathogen.AVR,
            "structure.str", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

pathogen.AVR.meta <- pathogen[match(rownames(pathogen.AVR), pathogen$`Isolate name`), 1:12]


genind <- adegenet::df2genind(pathogen.AVR, ploidy = 1, NA.char = "-9",
                              pop =  pathogen.AVR.meta$State)

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
pal = viridisLite::magma
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
if (length(pal) == 1) {
  PAL <- match.fun(pal)
  pal <- ggcompoplot::char2pal(post$population, PAL)
}

post$Decade <- floor(pathogen.AVR.meta$Year[match(post$sample, pathogen.AVR.meta$`Isolate name`)] / 10) * 10
# post$Decade <- factor(substr(post$Decade, 3,4),
#                       levels = c("70", "80", "90", "00", "10"))

grp2 <- max.col(dapc1$posterior)
names(grp2) <- rownames(dapc1$posterior)

xord <- dapc1$grp.coord |> dist() |> hclust()
xord <- xord$order
xord <- c(5,1,4,2,3)
grp2[] <- xord[grp2]

post$grp2 <-  grp2[as.character(post$sample)]

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
                 legend.position = "bottom",
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.3),
                                                           color = "black"),
                 panel.spacing.x = ggplot2::unit(c(0,0,0,0.2,
                                                   0,0.2,
                                                   0,0,0,0.2,
                                                   0,0,0,0,0.2,
                                                   0,0,0,0), "lines"),
                 axis.text.x.bottom = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggh4x::facet_nested(~Decade+grp2, 
                      space = "free_x",
                      scales = "free_x",) +
  ggplot2::scale_fill_manual(values = pal)

g1
