setwd("../Supplement")

source("~/scripts/r/utilities.r")

setup_data <- function(N, UPTO) {
  rice <-  readRDS(file = "rice.RDS") 
  pathogens <- readRDS(file = "pathogen.RDS") 
  
  #### this is just more manip
  # do some linear interpolations
  xat <- colnames(rice$Racreage.prop.known.bystate$AR$resistant)
  
  resistant <- lapply(rice$Racreage.prop.known.bystate,
                      \(x) x$resistant) |>
    lapply(\(z) apply(z, 1, \(r) approx(x = xat, 
                                        y = r, 
                                        xout = 1960:2018, rule = 2)$y) |>
             t()) |>
    lapply(\(z) {colnames(z) <- 1960:2018 ; return(z)})
  
  rm(xat)
  
  # use resistant "as-is" if we are not taking a rolling average
  if (N == 0 & UPTO == 0) {
  } else {
    resistant <- lapply(resistant,
                        \(x) {
                          x2 <- zoo::rollapply(t(x), list(seq(-N, -UPTO)), mean);
                          rownames(x2) <- colnames(x)[-c(1:N)];
                          x2 <- t(x2);
                          return(x2);
                        })
  }
  
  pathogens.right <- lapply(colnames(rice$Rgenes),
                            \(Rgene) purrr::map2_dbl(as.character(pathogens$Year), pathogens$State,
                                                     \(y, s) resistant[[s]][Rgene, y])) |>
    as.data.frame()
  
  
  colnames(pathogens.right) <- colnames(rice$Rgenes)
  
  pathogens[, -c(1:3)] <- lapply(pathogens[, -c(1:3)], as.factor)
  
  pathogens <- cbind(pathogens, pathogens.right)
  
  rm(pathogens.right)
  
  colnames(pathogens) <- stringr::str_replace_all(colnames(pathogens), pattern = " ", replacement = "_")
  colnames(pathogens) <- stringr::str_replace_all(colnames(pathogens), pattern = "-", replacement = "_")
  
  pathogens <- pathogens |>
    dplyr::select(-c(12:22))
  
  pathogens_short <- pathogens |>
    dplyr::mutate(dplyr::across(5:11, \(x) as.integer(x) - 1L)) |>
    dplyr::group_by(Year, State) |>
    dplyr::summarize(nTrials = dplyr::n(),
                     dplyr::across(dplyr::starts_with(c("ACE1", "AVR")), sum),
                     dplyr::across(dplyr::starts_with(c("Pi")), mean),
                     .groups = "drop")
  
  return(pathogens_short)
}

# replace values with the average of the preceeding N year
# excluding the current year
# test a parameter sweep
# N = num of years before current time
# UPTO = YEAR to stop at

lag_list <- build_timeseries_lag(N = 10, UPTO = 10)
parms <- lag_list$parms
mat <- lag_list$mat
rm(lag_list)

cognates <- list("AVR_Pib" = c("Pi_b"),
                 "AVR_Pizt" = c("Pi_z"),
                 "AVR_Pita1" = c("Pi_ta_Pi_ta2_Ptr"),
                 "AVR_Piks" = c("Pi_km", "Pi_ks"))

cognates.out <- cognates |>
  lapply(\(x) list())
## match it up with the AVR/R gene combinations and run all the models

for (cognate in names(cognates)) {
  outlist <- list()
  for (rown in 1:nrow(parms)) {
    N = parms$N[rown]
    UPTO = parms$UPTO[rown]
    Trial = parms$Trial[rown]
    
    data <- setup_data(N, UPTO)
    
    n_rgenes <- length(cognates[[cognate]])
    
    form <- paste0("cbind(", cognate, ", nTrials - ", cognate, ") ~ 1 + ",
                   "s(State, bs = 're') + s(Year, by = State, bs = 'cs') + ",
                   paste0(cognates[[cognate]], collapse = " + ")) |>
      as.formula()
    
    g <- mgcv::gam(form,
                   data = data,
                   family = "binomial", method = "REML")
  
    g.summary <- summary(g)
    
    range_ <- 2:(n_rgenes + 1)
    
    outlist[[rown]] <- list("coeff" = g.summary$p.coeff[range_] |> unname(),
                         "SE" = g.summary$se[range_] |> unname(),
                         "pval" = g.summary$p.pv[range_] |> unname(),
                         "AIC" = AIC(g),
                         "mod" = g,
                         "data" = data)
    names(outlist)[rown] <- Trial
  }
  cognates.out[[cognate]] <- outlist
}

# deal with the situation where there is more than one cognate R gene....
outlists.df <- lapply(cognates.out,
                      \(y) lapply(y, \(x) purrr::flatten_dbl(x[1:3]) |>
                                    matrix(ncol = 3, byrow = FALSE))) |>
  lapply(\(x) do.call(rbind, x)) |>
  lapply(as.data.frame) |>
  lapply(\(x) setNames(x, c("coeff", "SE", "pval"))) |>
  purrr::map2(cognates,
              \(df, cog) cbind(data.frame("Trial" = rep(parms$Trial, each = length(cog)),
                                          "R" = rep(cog, length(parms$Trial))),
                               df))

outlists.df.byRgene <- outlists.df |>
  lapply(\(x) dplyr::arrange(x, Trial)) |>
  lapply(\(x) split(x[, c("coeff", "SE", "pval")], x$R)) |>
  lapply(\(AVR) lapply(AVR, 
                       \(Rgene) lapply(Rgene, 
                                       \(x) {
                                        tmat <- mat |> t();
                                        tmat[upper.tri(tmat, diag = TRUE)] <- x;
                                        return(t(tmat));
                                        })))

outlists.df$AVR_Pita1 |> View()

for (AVR in names(cognates)) {
  for (Rgene in cognates[[AVR]]) {
    pheatmap::pheatmap(outlists.df.byRgene[[AVR]][[Rgene]][["coeff"]],
                       cellwidth = 30,
                       cellheight = 30,
                       fontsize_row = 8,
                       fontsize_col = 8,
                       main = paste("effect size of", Rgene, "on", AVR),
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       display_numbers = TRUE,
                       number_color = "black",
                       number_format = "%1.f",
                       # filename = "heatmap_pval.png",
                       # width = 10, height = 10
    )
    
    pheatmap::pheatmap(outlists.df.byRgene[[AVR]][[Rgene]][["pval"]] |> log10(),
                       cellwidth = 30,
                       cellheight = 30,
                       fontsize_row = 8,
                       fontsize_col = 8,
                       main =  paste("log10(pval) of", Rgene, "on", AVR),
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       display_numbers = TRUE,
                       number_color = "black",
                       number_format = "%.1f",
                       color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name =
                                                                           "RdYlBu"))(100)
                       # filename = "heatmap_pval.png",
                       # width = 10, height = 10
    )
    
  }
}







g <- cognates.out$AVR_Piks$`Trial_N-09_UPTO-05`$mod
data <- cognates.out$AVR_Piks$`Trial_N-09_UPTO-05`$data

# make combinations for the new data
nd <- expand.grid("Pi_km" = seq(0, 1, 0.05),
            "Pi_ks" = seq(0, 1, 0.1),      
            "State" = c("AR", "LA", "MS", "TX"),
            "Year" = c(1970, 1975,
                       1980, 1985,
                       1990, 1995,
                       2000, 2005,
                       2010, 2015)) |>
  dplyr::mutate(State = factor(State, levels = c("AR", "LA", "MS", "TX"))) |>
  dplyr::arrange(State)


nd <- expand.grid("Pi_km" = seq(0, 1, 0.05),
                  "Pi_ks" = seq(0, 1, 0.1),      
                  "State" = c("AR", "LA", "MS", "TX"),
                  "Year" = 1970:2015) |>
  dplyr::mutate(State = factor(State, levels = c("AR", "LA", "MS", "TX"))) |>
  dplyr::arrange(State)

x.terms <- predict(g, newdata = nd,
                   type = "iterms") 
colnames(x.terms) <- c("Rgene1", "Rgene2", "state", "AR", "LA", "MS", "TX")
x.terms <- cbind(matrix(rep(attr(x.terms, "constant"), nrow(x.terms)),
                        dimnames = list(c(),c("Intercept"))),
                 x.terms)

x.predicted <- cbind(nd, as.data.frame(x.terms)) |>
  dplyr::mutate(prob_AVR_Piks = plogis(rowSums(x.terms)))

ggplot2::ggplot(x.predicted,
                ggplot2::aes(x = Pi_km,
                             y = prob_AVR_Piks,
                             group = Year,
                             color = Year)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_color_viridis_c() +
  ggplot2::geom_line() +
  ggplot2::facet_grid(State~Pi_ks)

ggplot2::ggplot(x.predicted,
                ggplot2::aes(x = Pi_km,
                             y = prob_AVR_Piks,
                             color = Year,
                             group = Year)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_color_viridis_c() +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~State) +
  ggplot2::geom_point(data = data, 
                      ggplot2::aes(x = Pi_km,
                                   y = AVR_Piks/nTrials,
                                   size = nTrials),
                      alpha = 0.7,
                      inherit.aes = TRUE)



ggplot2::ggplot(x.predicted,
                ggplot2::aes(x = Year,
                             y = prob_AVR_Piks,
                             group = Pi_km,
                             color = Pi_km)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_color_viridis_c() +
  ggplot2::geom_line() +
  ggplot2::facet_grid(State~Pi_ks)







#####


g <- cognates.out$AVR_Pita1$`Trial_N-05_UPTO-05`$mod
data <- cognates.out$AVR_Pita1$`Trial_N-05_UPTO-05`$data

# make combinations for the new data
nd <- expand.grid("Pi_ta_Pi_ta2_Ptr" = seq(0, 1, 0.01),
                  "State" = c("AR", "LA", "MS", "TX"),
                  "Year" = c(1970, 1975,
                             1980, 1985,
                             1990, 1995,
                             2000, 2005,
                             2010, 2015)) |>
  dplyr::mutate(State = factor(State, levels = c("AR", "LA", "MS", "TX"))) |>
  dplyr::arrange(State)

x.terms <- predict(g, newdata = nd,
                   type = "iterms") 
colnames(x.terms) <- c("Rgene", "state", "AR", "LA", "MS", "TX")
x.terms <- cbind(matrix(rep(attr(x.terms, "constant"), nrow(x.terms)),
                        dimnames = list(c(),c("Intercept"))),
                 x.terms)

x.predicted <- cbind(nd, as.data.frame(x.terms)) |>
  dplyr::mutate(prob_AVR = plogis(rowSums(x.terms)))

ggplot2::ggplot(x.predicted,
                ggplot2::aes(x = Pi_ta_Pi_ta2_Ptr,
                             y = prob_AVR,
                             color = Year,
                             group = Year)) +
  ggplot2::scale_color_viridis_c() +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~State) +
  ggplot2::geom_point(data = data, 
                      ggplot2::aes(x = Pi_ta_Pi_ta2_Ptr,
                                   y = AVR_Pita1/nTrials,
                                   size = nTrials),
                      alpha = 0.7,
                      inherit.aes = TRUE)

