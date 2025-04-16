setwd("../../inputs")

build_timeseries_lag <- NULL

#' Build a time-series lag matrix for testing parameter combinations.
#'

#' @param N 
#' @param UPTO 
#'
#' @returns a list containing a `matrix` object of type NA_real_ with the lags
#' @export
#'
#' @examples
build_timeseries_lag <- function(N, UPTO) {
  mat <- matrix(NA_real_, nrow = N + 1, ncol = N + 1)
  
  parms <- which(lower.tri(mat,  diag = TRUE), arr.ind = TRUE) |>
    as.data.frame() |>
    dplyr::arrange(dplyr::desc(row), dplyr::desc(col)) |>
    setNames(c("N", "UPTO")) |>
    dplyr::mutate(N = N - 1L, UPTO = UPTO - 1L) |>
    dplyr::mutate("Trial" = paste0("Trial_", 
                                   "N-", stringr::str_pad(N, 2, pad = "0"),
                                   "_", 
                                   "UPTO-",stringr::str_pad(UPTO, 2, pad = "0")))
  
  rownames(mat) <- paste0("N-", stringr::str_pad(0:N, 2, pad = "0"))
  colnames(mat) <- paste0("UPTO-", stringr::str_pad(0:N, 2, pad = "0"))
  
  return(list(mat = mat, parms = parms))
}

R <- NULL
Pred <- NULL
proportion <- NULL
Year <- NULL
State <- NULL

setup_data <- function(N, UPTO) {
  rice  <- readRDS(file = "rice_smooth.RDS") |>
    tidyr::pivot_wider(names_from = R, values_from = Pred)
  pathogens <- readRDS(file = "pathogen.RDS") 
  
  resistant <- split(rice[, -c(1)], rice$State) |>
    lapply(\(x) {
      x2 <- data.matrix(x[, -c(1)]);
      rownames(x2) <- x$Year;
      return(t(x2));})
  
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
  
  # pivot from years in columns to R genes in columns
  resistant <- resistant |>
    lapply(data.frame, check.names = FALSE) |>
    lapply(tibble::rownames_to_column, "R") |>
    purrr::list_rbind(names_to = "State") |>
    tidyr::pivot_longer(-c(1:2), names_to = "Year", values_to = "proportion") |>
    tidyr::pivot_wider(names_from = R, values_from = proportion) |>
    dplyr::mutate(Year = as.double(Year))
  
  # join the host and pathogen data
  data <- dplyr::left_join(pathogens, resistant, by = c("State", "Year"),
                           relationship = "many-to-one")
  
  # fix the R gene names to work with glmmTMB
  colnames(data) <- stringr::str_replace_all(colnames(data), 
                                             pattern = " ", 
                                             replacement = "_")
  colnames(data) <- stringr::str_replace_all(colnames(data),
                                             pattern = "-",
                                             replacement = "_")
  
  # turn into binomial data
  pathogens_short <- data |>
    dplyr::mutate(dplyr::across(5:12, \(x) as.integer(x))) |>
    dplyr::group_by(Year, State,
                    dplyr::across(dplyr::starts_with(c("Pi", "pi")))) |>
    dplyr::summarize(nTrials = dplyr::n(),
                     dplyr::across(dplyr::starts_with(c("ACE1", "AVR")), sum),
                     dplyr::across(dplyr::starts_with(c("Pi")), mean),
                     .groups = "drop") |>
    dplyr::mutate(State = factor(State))
  
  return(pathogens_short)
}

cognates <- list("AVR_Pia" = c("Pi_a"),
                 "AVR_ii" = c("Pi_i"),
                 "ACE1" = c("pi_d"),
                 "AVR_Pib" = c("Pi_b"),
                 "AVR_Pizt" = c("Pi_z"),
                 "AVR_Pita1" = c("Pi_ta_Pi_ta2_Ptr"),
                 "AVR_Piks" = c("Pi_km", "Pi_ks"))

lag_list <- build_timeseries_lag(N = 14, UPTO = 14)
parms <- lag_list$parms
mat <- lag_list$mat
rm(lag_list)

cognates.out <- cognates |>
  lapply(\(x) list())
## match it up with the AVR/R gene combinations and run all the models

for (cognate in names(cognates)) {
  print(cognate)
  outlist <- list()
  for (rown in 1:nrow(parms)) {
    N = parms$N[rown]
    UPTO = parms$UPTO[rown]
    Trial = parms$Trial[rown]
    
    data <- setup_data(N, UPTO)
    
    n_rgenes <- length(cognates[[cognate]])
    
    form <- paste0("cbind(", cognate, ", nTrials - ", cognate, ") ~ 1 + ",
                   "(1 | State) + ",
                   paste0(cognates[[cognate]], collapse = " + ")) |>
      as.formula()
    
    
    g <- glmmTMB::glmmTMB(form,
                          data = data,
                          priors = data.frame(prior = c("normal(0, 5)"),
                                              class = "fixef",
                                              coef = ""),
                          family = "binomial",
                          REML = TRUE,
                          na.action = na.omit)

    g.summary <- summary(g)

    range_ <- 2:(n_rgenes + 1)
    
    outlist[[rown]] <- list("coeff" = g.summary$coefficients$cond[range_, 1] |> unname(),
                            "SE" =  g.summary$coefficients$cond[range_, 2] |> unname(),
                            "pval" =  g.summary$coefficients$cond[range_, 4] |> unname(),
                            "AIC" = AIC(g),
                            "mod" = g,
                            "data" = data |>
                              dplyr::mutate(Year_Start = Year - N,
                                            Year_End = Year - UPTO))

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

lapply(outlists.df, \(x) x$Trial[which.min(x$pval)])

outlists.best <- purrr::map2(cognates.out, lapply(outlists.df, \(x) x$Trial[which.min(x$pval)]),
                             \(x, y) x[[y]])

saveRDS(outlists.best, "../outputs/AVR-by-R_mods.RDS")

# lapply(outlists.best, \(x) performance::model_performance(x$mod))

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

newdata.list <- lapply(names(cognates),
                       \(x) list(list(Year = 1970),
                                 lapply(cognates[[x]],
                                        \(y) data.frame(x = seq(0, 1, 0.01)) |>
                                          setNames(y)) |>
                                   purrr::flatten())) |>
  lapply(purrr::flatten) |>
  lapply(expand.grid)

newdata.preds <- purrr::map2(outlists.best, newdata.list,
            \(x, df) predict(x$mod, newdata = df, re.form = ~ 0)) |>
  lapply(plogis) |>
  purrr::map2(newdata.list, \(x,y) cbind(y,x)) |>
  purrr::map2(names(cognates), \(x,y) {
    colnames(x)[ncol(x)] <- y;
    return(x);
  })

newdata.pred <- lapply(newdata.preds,
                       \(x) tidyr::pivot_longer(x, 2:(ncol(x) - 1),
                                                names_to = "R",
                                                values_to = "proportion") |>
                         tidyr::pivot_longer(2, names_to = "AVR", values_to = "prevalence")) |>
  purrr::list_rbind()


# get the metadata on the lags
meta.lag <- lapply(outlists.df, \(x) x$Trial[which.min(x$pval)])
meta.lag <- meta.lag[which(!names(meta.lag) %in% c("AVR_Piks", "AVR_Pib"))]
meta.lag <- lapply(meta.lag, 
                   \(x) paste0(substr(x, 9, 10), " to ",
                               substr(x, 17, 18), " years prior"))


g <- ggplot2::ggplot(newdata.pred |>
                  dplyr::filter(AVR != "AVR_Piks",
                                AVR != "AVR_Pib") |>
                  dplyr::arrange(AVR, R, proportion),
                ggplot2::aes(x = proportion,
                             y = prevalence)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~AVR,
                      labeller = ggplot2::labeller(AVR = function(x) {
                        x <- as.character(x);
                        return(sapply(x, \(z) paste0(z, " ~ ", cognates[[z]],
                                                     "\n", meta.lag[[z]])));
                      })) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 panel.spacing.x = ggplot2::unit(2, "lines")) +
  ggplot2::scale_y_continuous(name = "Predicted Proportion of Isolates with Cognate AVR Gene",
                              limits = c(0, 1)) +
  ggplot2::scale_x_continuous(name = "Proportion of Rice Acres with R Gene")

ggplot2::ggsave(g, filename = "../figures/Rgenes-AVR_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")


newdata.pred.piks <- newdata.preds$AVR_Piks |>
  dplyr::filter(Pi_ks %in% round(seq(0, 1, 1/6), 2))

meta.lag2 <- lapply(outlists.df, \(x) x$Trial[which.min(x$pval)])
meta.lag2 <- meta.lag2[which(names(meta.lag2) %in% c("AVR_Piks"))]
meta.lag2 <- lapply(meta.lag2, 
                   \(x) paste0(substr(x, 9, 10), " to ",
                               substr(x, 17, 18), " years prior")) |>
  unlist() |>
  unname()

g <- ggplot2::ggplot(newdata.pred.piks,
                ggplot2::aes(x = Pi_km,
                             y = AVR_Piks)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~Pi_ks, labeller = ggplot2::labeller(Pi_ks = function(x) {
    x <- as.character(x)
    return(paste0("AVR_Piks ~ Pi_km | Pi-ks: ", x, "\n", meta.lag2))
  })) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14, family = "Arial", color = "black"),
                 legend.position = "inside",
                 legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                 panel.spacing.x = ggplot2::unit(2, "lines")) +
  ggplot2::scale_y_continuous(name = "Predicted Proportion of Isolates with Cognate AVR Gene",
                              limits = c(0, 1)) +
  ggplot2::scale_x_continuous(name = "Proportion of Rice Acres with R Gene")

ggplot2::ggsave(g, filename = "../figures/Rgenes-AVR-Piks_over_time.png",
                width = 6,
                height = 3,
                dpi = 600,
                scale = 2,
                units = "in")
