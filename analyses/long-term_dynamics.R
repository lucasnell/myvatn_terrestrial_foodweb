

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(grid)
    library(pbmcapply)
})


dir <- sprintf("~/Box Sync/Iceland Food Web Model/Results/Figures_%s/", Sys.Date())

if (!dir.exists(dir)) dir.create(dir)

#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'




full_pulse_df <- read_csv("~/Box Sync/Iceland Food Web Model/Results/sim_combinations.csv.gz",
                         col_types = cols(
                             w = "f",
                             pool = "f",
                             .default = "d")) %>%
        mutate_at(vars(b, f), ~ factor(signif(., digits = 6)))

make_vec <- function(cname) rep(signif(par_combs[[cname]], digits = 6), each = 6)

par_combs <- expand.grid(w = seq(10, 25, length.out = 25),
                         b = seq(0.1, 40, length.out = 25),
                         f = seq(8e-3, 8, length.out = 25),
                         mM = par_estimates$mM[1] * c(0.5, 1, 2),
                         hM = par_estimates$hM[1] * c(0.5, 1, 2),
                         # Plant / herbivore uptake rates:
                         aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                         aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    mutate_at(vars(b, f, mM, hM, aDV, aPH), ~ signif(., digits = 6))


pulse_df <- full_pulse_df %>%
    mutate(mM = make_vec("mM"),
           hM = make_vec("hM"),
           aDV = make_vec("aDV"),
           aPH = make_vec("aPH")) %>%
    filter(mM == signif(par_estimates$mM[1], digits = 6),
           hM == signif(par_estimates$hM[1], digits = 6),
           aDV == signif(par_estimates$aDV[1], digits = 6),
           aPH == signif(par_estimates$aPH[1], digits = 6)) %>%
    select(-mM, -hM, -aDV, -aPH)



heat_plot <- function(.col, .facet_lab = NULL,
                      .pool = "detritivore", .df = NULL, .w = NULL, .title = NULL,
                      .legend_title = NULL,
                      .legend_breaks = waiver(), .legend_labels = waiver(),
                      flip_sign = FALSE, flip_color = FALSE,
                      keep_x = TRUE, keep_y = TRUE, .pal_opt = "C") {
    if (is.null(.df)) .df <- pulse_df
    if (is.null(.w)) .w <- .df$w %>% paste() %>% as.numeric() %>% median()
    .col <- enquo(.col)
    # if (is.null(.title)) .title <- .pool
    .df <- .df %>%
        filter(pool == .pool, w == .w) %>%
        mutate(f = f %>% paste() %>% as.numeric())
    if (flip_sign) {
        .df <- .df %>% mutate_at(vars(!!.col), ~ . * -1)
    }
    p <- .df %>%
        ggplot(aes(f, area)) +
        geom_tile(aes(fill = !!.col)) +
        geom_contour(aes(z = !!.col), color = "white", bins = 7) +
        scale_fill_viridis_c(.legend_title, option = .pal_opt, breaks = .legend_breaks,
                             labels = .legend_labels) +
        # ggtitle(.title) +
        ylab(expression("Total midge input (" * g ~ N ~ m^{-2} * ")")) +
        scale_x_continuous("Midge accessibility", breaks = c(0, 4, 8)) +
        theme(legend.title = element_text(size = 10),
              legend.key.width = grid::unit(0.02, "npc"),
              plot.title = element_text(hjust = 0, size = 12),
              axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12),
              plot.margin = margin(0,0,0,0),
              legend.margin = margin(0,0,0,0),
              legend.text = element_text(size = 10))
    if (!is.null(.facet_lab)) {
        p <- p  +
            geom_text(data = tibble(f = 0.008, area = 750, lab = .facet_lab),
                      aes(label = lab), hjust = 0, vjust = 0, size = 12 / 2.835)
    }
    if (!keep_x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    if (!keep_y) {
        p <- p + theme(axis.title.y = element_blank())
    }
    return(p)
}



hm1 <- heat_plot(cum_pos_loss_V, .pal_opt = "B", .facet_lab = "a",
                 .legend_breaks = c(0.1, 1.5),
                 keep_y = FALSE, keep_x = FALSE,
                 .legend_title = "Top-down\nintensification")
hm2 <- heat_plot(cum_neg_loss_V, .pal_opt = "D", flip_sign = TRUE, .facet_lab = "b",
                 .legend_breaks = c(0.1, 0.7),
                 .legend_title = "Top-down\nalleviation", keep_x = FALSE)
hm3 <- heat_plot(cum_gain_V, .legend_title = "Bottom-up", .facet_lab = "c",
                 keep_y = FALSE,
                 .pal_opt = "C") +
    scale_fill_gradient("Bottom-up", high = "deepskyblue", low = "black",
                        breaks = c(1, 5))
# heat_plot(cum_gain_H, "herbivore", .title = "Enhancement of bottom-up effect", flip_color = TRUE)





hmg1 <- ggplotGrob(hm1)
hmg2 <- ggplotGrob(hm2)
hmg3 <- ggplotGrob(hm3)

hm <- rbind(hmg1, hmg2, hmg3, size = "first")
hm$widths <- unit.pmax(hmg1$widths, hmg2$widths, hmg3$widths)



# grid.newpage(); grid.draw(hm)


pdf(file = paste0(dir, "5-heatmaps.pdf"), width = 4.5, height = 9)
grid.newpage()
grid.draw(hm)
dev.off()








pulse_df_fig6 <- full_pulse_df %>%
    mutate(mM = make_vec("mM"),
           hM = make_vec("hM"),
           aDV = make_vec("aDV"),
           aPH = make_vec("aPH")) %>%
    filter(aDV == signif(par_estimates$aDV[1], digits = 6),
           aPH == signif(par_estimates$aPH[1], digits = 6),
           w == 20,
           pool == "soil") %>%
    select(-aDV, -aPH, -w, -pool, -b) %>%
    mutate(mM = factor(mM, levels = sort(unique(mM)),
                       labels = paste(c("low", "mid", "high"), "midge\ndecay rate")),
           hM = factor(hM, levels = sort(unique(hM)),
                       labels = paste(c("low", "mid", "high"),
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
               grepl("^cum_loss_", variable) ~ "TD total",
               TRUE ~ "BU total"
           ), levels = c("BU total", "TD total",
                         "TD intensification", "TD alleviation")),
           f = as.numeric(paste(f))) %>%
    select(-variable) %>%
    filter(mM == "mid midge\ndecay rate") %>%
    filter(f %in% quantile(f, c(0.1, 0.5, 0.9))) %>%
    mutate(f = factor(f, levels = sort(unique(f)),
                      labels = sprintf("%s midge\naccessibility",
                                       c("low", "mid", "high")))) %>%
    mutate(id = interaction(direction, type)) %>%
    filter(area < 350)



td_bu_avail_plot <- pulse_df_fig6 %>%
    filter(pool == "detritivore") %>%
    ggplot(aes(area, value)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_line(aes(color = direction, linetype = type, group = id), size = 1) +
    facet_grid(hM ~ f) +
    scale_color_manual(NULL, values = c(color_pal()[1], "gray60"), guide = FALSE) +
    scale_linetype_manual(NULL, values = c(1, 1, 3:2)) +
    guides(linetype = guide_legend(keywidth = 2, nrow = 2, keyheight = 0.6,
                                   override.aes = list(color = c(color_pal()[1],
                                                                 rep("gray60", 3))),
                                   byrow = TRUE)) +
    xlab(expression("Total midge input (" * g ~ N ~ m^{-2} * ")")) +
    ylab("Cumulative midge effect") +
    theme(strip.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.box = "vertical",
          legend.title = element_text(size = 11, margin = margin(0,0,0,b=-4)),
          legend.text = element_text(size = 10)) +
    NULL


ggsave(filename = paste0(dir, "6-td_bu_avail.pdf"), td_bu_avail_plot, width = 7, height = 5)






# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================



parlist <- par_estimates %>%
    filter(V==1, H==1, R==1, iN == 10) %>%
    as.list()
# Top-down
V_gain <- function(V, D, aDV) {
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
# Bottom-up
V_loss <- function(V, R, H, M, f, hM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    ((aR*V*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / V
}
H_gain <- function(P, H, aPH) {
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}
H_loss <- function(H, R, V, M, f, hM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    ((aR*H*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / H
}




# Runs simulations for one combination of parameters for Fig. 7

one_combo_fig7 <- function(row_i) {
    .w <- 20
    .b <- row_i$b
    .other_pars <- as.list(unlist(row_i))
    .other_pars$b <- NULL

    fw <- food_web(tmax = 250, s = 10, b = .b, w = .w, other_pars = .other_pars)

    fw <- fw %>%
        spread(pool, N) %>%
        mutate(Vg = V_gain(detritivore, detritus, par_estimates$aDV[1]),
               Vl = V_loss(detritivore, predator, herbivore, midge, .other_pars$f,
                           .other_pars$hM),
               Hg = H_gain(plant, herbivore, par_estimates$aPH[1]),
               Hl = H_loss(herbivore, predator, detritivore, midge, .other_pars$f,
                           .other_pars$hM)) %>%
        gather("pool", "N", soil:midge) %>%
        filter(pool != "midge") %>%
        mutate(pool = factor(pool,
                             levels = c("soil", "detritus", "plant",
                                        "detritivore", "herbivore", "predator"))) %>%
        arrange(pool, time) %>%
        group_by(pool) %>%
        summarize(cum_pos_loss_V = sum(Vl[Vl > Vl[1]] - Vl[1]),
                  cum_neg_loss_V = sum(Vl[Vl < Vl[1]] - Vl[1]),
                  cum_pos_loss_H = sum(Hl[Hl > Hl[1]] - Hl[1]),
                  cum_neg_loss_H = sum(Hl[Hl < Hl[1]] - Hl[1]),
                  cum_gain_V = sum(Vg[Vg > Vg[1]] - Vg[1]),
                  cum_gain_H = sum(Hg[Hg > Hg[1]] - Hg[1])) %>%
        ungroup() %>%
        mutate(b = .b, mM = .other_pars$mM, hM = .other_pars$hM) %>%
        select(b, everything())
    return(fw)
}

pulse_df_fig7 <- crossing(b = seq(5, 50, length.out = 10),
                      f = 3,
                      mM = par_estimates$mM[1] * c(0.5, 1, 2),
                      hM = par_estimates$hM[1] * c(0.5, 1, 2)) %>%
    split(row(.)[,1]) %>%
    pbmclapply(one_combo_fig7, mc.cores = parallel::detectCores()) %>%
    bind_rows() %>%
    filter(pool %in% c("detritivore", "herbivore")) %>%
    mutate(mM = factor(mM, levels = sort(unique(mM)),
                       labels = paste(c("low", "mid", "high"), "midge\ndecay rate")),
           hM = factor(hM, levels = sort(unique(hM)),
                       labels = paste(c("low", "mid", "high"),
                                      "midge\nhandling time"))) %>%
    group_by(b, pool, mM, hM) %>%
    summarize(cum_pos_loss = ifelse(pool[1] == "detritivore", cum_pos_loss_V[1],
                                    cum_pos_loss_H[1]),
              cum_neg_loss = ifelse(pool[1] == "detritivore", cum_neg_loss_V[1],
                                    cum_neg_loss_H[1]),
              cum_loss = cum_pos_loss + cum_neg_loss,
              cum_gain = ifelse(pool[1] == "detritivore", cum_gain_V[1],
                                cum_gain_H[1])) %>%
    ungroup() %>%
    arrange(b)



td_bu_plot <- pulse_df_fig7 %>%
    mutate(b = b * 20) %>%
    rename(area = b) %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_path(aes(color = pool), size = 1) +
    geom_point(aes(size = area, color = pool)) +
    geom_text(data = tibble(cum_gain = c(5, 1),
                            cum_loss = rep(0.5, 2),
                            pool = sort(unique(pulse_df_fig7$pool)),
                            hM = sort(unique(pulse_df_fig7$hM))[3],
                            mM = sort(unique(pulse_df_fig7$mM))[1]),
              aes(label = pool, color = pool), hjust = 0, vjust = 1, size = 10 / 2.835) +
    facet_grid(mM ~ hM) +
    scale_color_manual(values = color_pal()[1:2], guide = FALSE) +
    scale_size_continuous(range = c(0.5, 4),
                          breaks = c(100, 500, 1000)) +
    guides(size = guide_legend(expression("Total midge input (" * g ~ N ~ m^{-2} * ")"),
                               title.position = "top")) +
    xlab("Cumulative midge effect on BU") +
    ylab("Cumulative midge effect on TD") +
    theme(strip.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.title = element_text(size = 11, margin = margin(0,0,0,b=-4)),
          legend.text = element_text(size = 10)) +
    NULL



# ggsave(filename = paste0(dir, "7-top_vs_bottom.pdf"), td_bu_plot, width = 7, height = 5)

