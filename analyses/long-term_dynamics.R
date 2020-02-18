

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(grid)
    # library(pbmcapply)
})


dir <- sprintf("~/Box Sync/Iceland Food Web Model/Results/Figures_%s/", Sys.Date())

if (!dir.exists(dir)) dir.create(dir)

#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'





pulse_df_fig6 <- read_csv(paste0("~/Box Sync/Iceland Food Web Model/Results/",
                                 "sim_combinations.csv.gz"),
                          col_types = cols(pool = "f", .default = "d")) %>%
    filter(aDV == par_estimates$aDV[1],
           aPH == par_estimates$aPH[1],
           pool == "soil",
           mM == median(mM),
           f %in% range(f),
           hM %in% range(hM)) %>%
    mutate(area = w * b,
           f = factor(f, levels = sort(unique(f)),
                      labels = sprintf("%s midge\nexploitation",
                                       c("low", "high"))),
           hM = factor(hM, levels = sort(unique(hM)),
                       labels = paste(c("low", "high"),
                                      "midge\nhandling time"))) %>%
    select(area, f, mM, hM,
           cum_pos_loss_V, cum_pos_loss_H,
           cum_neg_loss_V, cum_neg_loss_H,
           cum_gain_V, cum_gain_H) %>%
    mutate(cum_loss_V = cum_pos_loss_V + cum_neg_loss_V,
           cum_loss_H = cum_pos_loss_H + cum_neg_loss_H,
           # Changing to magnitudes
           cum_neg_loss_V = abs(cum_neg_loss_V),
           cum_neg_loss_H = abs(cum_neg_loss_H)) %>%
    gather("variable", "value", cum_pos_loss_V:cum_loss_H) %>%
    mutate(pool = factor(ifelse(grepl("V$", variable), "detritivore", "herbivore")),
           direction = factor(ifelse(grepl("loss", variable), "top-down", "bottom-up")),
           type = factor(case_when(
               grepl("^cum_pos_loss", variable) ~ "TD intensification",
               grepl("^cum_neg_loss", variable) ~ "TD alleviation",
               grepl("^cum_loss_", variable) ~ "TD net",
               TRUE ~ "BU total"
           ), levels = c("BU total", "TD net",
                         "TD intensification", "TD alleviation"))) %>%
    select(-variable) %>%
    mutate(id = interaction(direction, type)) %>%
    filter(area < 500) %>%
    identity()



fig6_labs <- list(bquote(italic('BU'['total'])),
                  bquote(italic('TD'['net'])),
                  bquote(italic('TD'['intensification'])),
                  bquote(italic('TD'['alleviation'])))



td_bu_avail_plot <- pulse_df_fig6 %>%
    filter(pool == "detritivore") %>%
    ggplot(aes(area, value)) +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_line(aes(color = direction, linetype = type, group = id), size = 1) +
    geom_text(data = tibble(area =  0, value = max(pulse_df_fig6$value),
                            hM = pulse_df_fig6$hM %>% unique() %>%
                                sort() %>% rep(each = 2),
                            f = pulse_df_fig6$f %>% unique() %>%
                                sort() %>% rep(2),
                            labs = letters[1:4]),
              aes(label = labs), hjust = 0, vjust = 1,
              size = 12 / 2.835) +
    facet_grid(hM ~ f) +
    scale_color_manual(NULL, values = c(color_pal()[1], "gray60"), guide = FALSE) +
    scale_linetype_manual(NULL, values = c(1, 1, 3:2),
                          breaks = c("BU total", "TD net",
                                     "TD intensification", "TD alleviation"),
                          labels = fig6_labs) +
    guides(linetype = guide_legend(keywidth = 2, nrow = 2, keyheight = 0.6,
                                   override.aes = list(color = c(color_pal()[1],
                                                                 rep("gray60", 3))),
                                   byrow = TRUE)) +
    xlab(expression("Total midge input (" * g ~ N ~ m^{-2} * ")")) +
    scale_y_continuous("Total midge effect", limits = c(-0.4, 4.3)) +
    theme(strip.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.box = "vertical",
          legend.title = element_text(size = 11, margin = margin(0,0,0,b=-4)),
          legend.text = element_text(size = 10),
          strip.text.y = element_text(vjust = 1)) +
    NULL



cairo_pdf(filename = paste0(dir, "5-td_bu_avail.pdf"), width = 5, height = 5)
td_bu_avail_plot
dev.off()


