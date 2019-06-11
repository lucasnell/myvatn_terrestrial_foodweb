

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
                         f = seq(8e-3, 8e-2, length.out = 25),
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
        ylab("Total midge input") +
        scale_x_continuous("Midge accessibility", breaks = c(0.04, 0.08)) +
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





hm1 <- heat_plot(cum_pos_loss_V, .pal_opt = "B", #.facet_lab = "a",
                 .legend_breaks = c(0.1, 1.5),
                 # keep_y = FALSE, keep_x = FALSE,
                 .legend_title = "Top-down\nintensification")
hm2 <- heat_plot(cum_neg_loss_V, .pal_opt = "D", flip_sign = TRUE, # .facet_lab = "b",
                 .legend_breaks = c(0.01, 0.075),
                 .legend_title = "Top-down\nalleviation") # , keep_x = FALSE)
hm3 <- heat_plot(cum_gain_V, .legend_title = "Bottom-up",  #.facet_lab = "c",
                 #keep_y = FALSE,
                 .pal_opt = "C") +
    scale_fill_gradient("Bottom-up", high = "deepskyblue", low = "black",
                        breaks = c(1, 5))
# heat_plot(cum_gain_H, "herbivore", .title = "Enhancement of bottom-up effect", flip_color = TRUE)


library(cowplot)

plot_grid(hm1, hm2, hm3, trans_p4, align = "hv", axis = "lrb")


hmg1 <- ggplotGrob(hm1)
hmg2 <- ggplotGrob(hm2)
hmg3 <- ggplotGrob(hm3)

hm <- rbind(hmg1, hmg2, hmg3, size = "first")
hm$widths <- unit.pmax(hmg1$widths, hmg2$widths, hmg3$widths)



# grid.newpage()
# grid.draw(hm)


# pdf(file = paste0(dir, "5-heatmaps.pdf"), width = 4.5, height = 9)
# grid.newpage()
# grid.draw(hm)
# dev.off()

pdf(file = paste0(dir, "5-heatmaps2.pdf"), width = 7, height = 5)
plot_grid(hm1, hm2, hm3, trans_p4, align = "hv", axis = "lrb", labels = letters[1:4])
dev.off()



parlist <- par_estimates %>%
    filter(V==1, H==1, R==1, iN == 10) %>%
    as.list()
V_gain <- function(V, D, aDV) {
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
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

# Runs simulations for one combination of parameters
one_combo <- function(row_i) {
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




pulse_df2 <- crossing(b = seq(5, 50, length.out = 10),
                      f = 3,
                      mM = par_estimates$mM[1] * c(0.5, 1, 2),
                      hM = par_estimates$hM[1] * c(0.5, 1, 2)) %>%
    split(row(.)[,1]) %>%
    pbmclapply(one_combo, mc.cores = parallel::detectCores()) %>%
    bind_rows() %>%
    filter(pool %in% c("detritivore", "herbivore")) %>%
    mutate(mM = factor(mM, levels = sort(unique(mM)),
                       labels = paste(c("low", "mid", "high"), "midge\ndecay rate")),
           hM = factor(hM, levels = sort(unique(hM)),
                       labels = paste(c("high", "mid", "low"),
                                      "midge\naccessibility"))) %>%
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



td_bu_plot33 <- pulse_df2 %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_path(aes(color = pool), size = 1) +
    geom_point(aes(size = b, color = pool)) +
    geom_text(data = tibble(cum_gain = c(5, 1),
                            cum_loss = rep(0.3, 2),
                            pool = sort(unique(pulse_df2$pool)),
                            hM = sort(unique(pulse_df2$hM))[3],
                            mM = sort(unique(pulse_df2$mM))[1]),
              aes(label = pool, color = pool), hjust = 0, vjust = 1, size = 10 / 2.835) +
    facet_grid(mM ~ hM) +
    scale_color_manual(values = color_pal()[1:2], guide = FALSE) +
    scale_size_continuous("Midge pulse\nintensity", range = c(0.5, 4),
                          breaks = c(1, 5, 25)) +
    xlab("Bottom-up effect") +
    ylab("Total top-down effect") +
    # theme(legend.position = c(0.85, 0.25), legend.background = element_blank()) +
    NULL



# ggsave(filename = paste0(dir, "6-top_vs_bottom_3x3.pdf"), td_bu_plot33, width = 8, height = 5)



td_bu_plot22 <- pulse_df2 %>%
    filter(!grepl("^mid", mM), !grepl("^mid", hM)) %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_path(aes(color = pool), size = 1) +
    geom_point(aes(size = b, color = pool)) +
    geom_text(data = tibble(cum_gain = c(5, 1),
                            cum_loss = rep(0.2, 2),
                            pool = sort(unique(pulse_df2$pool)),
                            hM = sort(unique(pulse_df2$hM))[3],
                            mM = sort(unique(pulse_df2$mM))[1]),
              aes(label = pool, color = pool), hjust = 0, vjust = 1, size = 10 / 2.835) +
    facet_grid(mM ~ hM) +
    scale_color_manual(values = color_pal()[1:2], guide = FALSE) +
    scale_size_continuous("Midge pulse\nintensity", range = c(0.5, 4),
                          breaks = c(1, 5, 25)) +
    xlab("Bottom-up effect") +
    ylab("Total top-down effect") +
    # theme(legend.position = c(0.85, 0.25), legend.background = element_blank()) +
    NULL


# ggsave(filename = paste0(dir, "6-top_vs_bottom_2x2.pdf"), td_bu_plot22, width = 8, height = 5)
