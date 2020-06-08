# ======================================*
# Plots for transient dynamics
# ======================================*



#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'




# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(viridisLite)
    library(grid)
    library(egg)
})


# This sets plotting device on LAN computer:
if (Sys.getenv("RSTUDIO") == "1" && file.exists(".Rprofile")) {
    source(".Rprofile")
}






dir <- sprintf("~/Box Sync/Iceland Food Web Model/Results/Figures_%s/", Sys.Date())

if (!dir.exists(dir)) dir.create(dir)




# ================================================================================*
# ================================================================================*

# Figure 2 ----

# ================================================================================*
# ================================================================================*



fig2_caption <- paste("Time series of proportional changes in N content among",
                      "(a) upper trophic level (herbivore, detritivore, predator) and",
                      "(b) lower trophic level (detritus, soil, plant) pools.",
                      "Lines indicate the N content at time $t$ relative to that at",
                      "time 0 [i.e. $(N_t − N_0) / N_0$].",
                      "The lower black bar represents the duration of the midge N pulse.",
                      "Parameter values are as in Table 1, with predator exploitation",
                      "$q = 3$, pulse duration $w = 20$, and pulse rate $b = 50$.")


upper_levels <- c("detritivore", "herbivore", "predator")



fig2_sim <- food_web(tmax = 100, s = 10, b = 50, w = 20,
                     other_pars = list(q = 3)) %>%
    filter(pool != "midge") %>%
    mutate(pool = droplevels(pool),
           level = ifelse(pool %in% upper_levels, 0, 1) %>%
               factor(levels = 0:1,
                      labels = paste(c("Upper", "Lower"), "trophic levels"))) %>%
    group_by(pool) %>%
    mutate(N = (N - N[1]) / (N[1])) %>%
    ungroup()



fig2_p <- fig2_sim %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_segment(data = tibble(time = 10, time2 = 10+20, N = -0.25),
                 aes(xend = time2, yend = N), size = 1.5) +
    geom_line(aes(color = pool), size = 1) +
    scale_y_continuous("Proportional change in N", breaks = c(0, 2, 4)) +
    scale_x_continuous("Time (days)") +
    facet_wrap(~ level, nrow = 2) +
    scale_color_manual(values = color_pal()[c(4:6, 1:3)]) +
    theme(legend.position = "none",
          panel.spacing = unit(1.5, "lines"),
          axis.title.y = element_text(margin = margin(0,0,0,r=12)),
          axis.title.x = element_text(margin = margin(0,0,0,t=12)),
          strip.text = element_text(size = 12, margin = margin(0,0,0,b=10))) +
    geom_text(data = tibble(
        pool = sort(unique(fig2_sim$pool[fig2_sim$pool != "midge"])),
        time =  c( 43,  15,  30,
                   30,  31.5,  48),
        N =     c(2.1, 2.3, 0.8,
                  1.9, 0.25, 3),
        level = factor(rep(1:0, each = 3), levels = 0:1,
                       labels = paste(c("Upper", "Lower"), "trophic levels"))),
        aes(label = pool, color = pool), size = 10 / 2.835) +
    geom_text(data = tibble(time = 20, N = -0.35,
                            level = factor(0, levels = 0:1,
                                           labels = paste(c("Upper", "Lower"),
                                                          "trophic levels"))),
              label = "pulse", size = 10 / 2.835, hjust = 0.5, vjust = 1,
              color = "black") +
    geom_text(data = tibble(time =  rep(-5, 2),
                            N = rep(5, 2),
                            level = factor(paste(c("Upper", "Lower"), "trophic levels")),
                            labs = sprintf("(%s)", letters[1:2])),
              aes(label = labs), hjust = 1, vjust = 0, size = 14 / 2.835) +
    coord_cartesian(ylim = c(-0.5, max(fig2_sim$N)),
                    xlim = c(0, 100),
                    clip = "off")


cairo_pdf(file = paste0(dir, "2-N_timeseries.pdf"), width = 5, height = 5)
fig2_p
dev.off()






# ================================================================================*
# ================================================================================*

# Figure 3 ----

# ================================================================================*
# ================================================================================*


fig3_caption <- paste("Time series of proportional changes in N content for",
                      "detritivore and herbivore pools (indicated by the left",
                      "y-axis).",
                      "This is shown for when",
                      "(a) midges only go to the detritus pool,",
                      "(b) midges go to both detritus and predator pools, and",
                      "(c) midges only go to the predator pool.",
                      "The blue shaded regions represent the proportional change in",
                      "N for predators, which is indicated by the right y-axis.",
                      "Different y-axis scales are used to aid visualization of",
                      "the transient dynamics of the herbivores and detritivores,",
                      "which had a weaker response to the pulse than predators for",
                      "the selected parameter values.",
                      "The lower black bar represents the duration of the midge",
                      "N pulse.",
                      "Parameter values are as in Table 1, with",
                      "predator exploitation $q = 3$,",
                      "pulse duration $w = 20$, and pulse rate $b = 20$.")




do_sim_fig3 <- function(.midges_not_to) {

    # .midges_not_to = "none"
    # rm(.midges_not_to)

    stopifnot(length(.midges_not_to) == 1 && .midges_not_to %in% c("none", "D", "X"))

    sim <- food_web(tmax = 100, s = 10, b = 20, w = 20,
                                other_pars = list(q = 3),
                    midges_not_to = .midges_not_to) %>%
        filter(pool %in% c(upper_levels, "midge")) %>%
        mutate(pool = droplevels(pool)) %>%
        mutate(pool = factor(paste(pool), levels = c(upper_levels, "midge")),
               not_to = .midges_not_to)

    return(sim)
}


fig3_df <- tibble(.midges_not_to = c("X", "none", "D")) %>%
    pmap_dfr(do_sim_fig3) %>%
    group_by(not_to, pool) %>%
    mutate(N_rel = (N - N[1]) / N[1]) %>%
    mutate(N_rel = ifelse(is.nan(N_rel), 0, N_rel)) %>%
    ungroup()



fig3_panel <- function(.midges_not_to,
                       .mult = 3.5,
                       .ylims = c(-0.303, 1.25)) {

    # .midges_not_to = "X"; .mult = 3.5; .ylims = c(-0.3, 1.25)
    # rm(.midges_not_to)

    stopifnot(length(.midges_not_to) == 1 && .midges_not_to %in% c("none", "D", "X"))

    .title <- case_when(.midges_not_to == "none" ~ "Midges to detritus and predators",
                        .midges_not_to == "D" ~ "Midges to predators only",
                        .midges_not_to == "X" ~ "Midges to detritus only",
                        TRUE ~ "ERROR")

    mt_combo <- fig3_df %>%
        filter(not_to == .midges_not_to) %>%
        select(-not_to)



    mt_plot <- mt_combo %>%
        filter(!pool %in% c("midge", "predator")) %>%
        mutate(pool = droplevels(pool)) %>%
        ggplot(aes(time, N_rel)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_segment(data = tibble(time = 10, time2 = 10+20, N_rel = -0.15),
                     aes(xend = time2, yend = N_rel), size = 1.5) +
        geom_ribbon(data = mt_combo %>%
                        filter(pool == "predator") %>%
                        mutate(N_rel = N_rel / .mult),
                    aes(ymin = 0, ymax = N_rel), fill = color_pal(0.5)[3]) +
        geom_line(aes(color = pool), size = 0.75) +
        scale_x_continuous("Time (days)") +
        scale_color_manual(NULL, values = color_pal()) +
        ggtitle(.title) +
        theme(legend.position = "none",
              strip.text.y = element_text(face = "plain", size = 10, angle = 270,
                                          margin = margin(l = 4)),
              strip.text.x = element_text(face = "plain", size = 10,
                                          margin = margin(b = 4)),
              panel.spacing.y = unit(0, "lines"),
              panel.spacing.x = unit(1.5, "lines"),
              plot.margin = margin(0,0,t=12,b=12),
              plot.title = element_text(size = 12, hjust = 0.5, face = "plain",
                                        margin = margin(0,0,0,b = 10))) +
        scale_y_continuous("Proportional change in N", limits = .ylims,
                           sec.axis = sec_axis(~ . * .mult,
                                               breaks = c(0, 2, 4),
                                               name = "Proportional change in N (predator)"),
                           breaks = c(0, 0.5, 1)) +
        NULL

    if (.midges_not_to == "none") {

        mt_plot <- mt_plot +
            geom_text(data = tibble(time =  65,
                                    N_rel =     1.2),
                      label = "predator", color = color_pal()[3], hjust = 1, vjust = 1,
                      size = 10 / 2.835) +
            geom_text(data = tibble(pool = factor(c("detritivore", "herbivore")),
                                    time =  c(92, 100),
                                    N_rel =     c(0.35, -0.15)),
                      aes(label = pool, color = pool), hjust = 1, size = 10 / 2.835) +
            geom_text(data = tibble(time = 20, N_rel = -0.2),
                      label = "pulse", size = 10 / 2.835, hjust = 0.5, vjust = 1,
                      color = "black") +
            theme(axis.title.y.left = element_text(margin = margin(0,0,0,r=12),
                                                   size = 11),
                  axis.title.y.right = element_text(margin = margin(0,0,0,l=12),
                                                    size = 11))
    } else {

        mt_plot <- mt_plot +
            theme(axis.title.y = element_blank())

    }

    if (.midges_not_to != "D") {
        mt_plot <- mt_plot +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    } else {
        mt_plot <- mt_plot +
            theme(axis.title.x = element_text(margin = margin(0,0,0,t=4), size = 11))
    }

    return(mt_plot)

}



fig3_panel_list <- tibble(.midges_not_to = c("X", "none", "D")) %>%
    pmap(fig3_panel) %>%
    map(~ .x + theme(axis.line.y.left = element_line(size = 1)))


cairo_pdf(file = paste0(dir, "3-N_midge_attack.pdf"), width = 3.5, height = 6.5)
ggarrange(plots = fig3_panel_list, ncol = 1, labels = sprintf("(%s)", letters[1:3]),
          label.args = list(gp = gpar(fontsize = 14, fontface =  "plain"),
                            vjust = 2, hjust = -1.5))
dev.off()





# ================================================================================*
# ================================================================================*

# Figure 4 ----

# ================================================================================*
# ================================================================================*


fig4_caption <- paste("Time series of the top-down (\"TD\") and bottom-up (\"BU\")",
                      "effects on detritivore and herbivore pools.",
                      "This is shown for when",
                      "(a) midges only go to the detritus pool,",
                      "(b) midges go to both detritus and predator pools, and",
                      "(c) midges only go to the predator pool.",
                      "The gray line indicates the top-down effects on both",
                      "detritivore and herbivore pools, as they are identical",
                      "in the model.",
                      "Parameter values are as in Table 1, with",
                      "predator exploitation $q = 3$,",
                      "pulse duration $w = 20$, and pulse rate $b = 20$.")




parlist <- par_estimates %>%
    filter(V==1, H==1, X==1, iI == 10) %>%
    as.list()

V_gain <- function(V, D, aDV = parlist[["aDV"]]) {
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}

V_loss <- function(V, X, H, M, q, hM = parlist[["hM"]], .no_M_to_X = FALSE) {
    aX <- parlist[["aX"]]
    hX <- parlist[["hX"]]
    if (!.no_M_to_X) {
        Vl <- ((aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / V
    } else {
        Vl <- ((aX*V*X)/(1 + aX*hX*(V + H))) / V
    }
    return(Vl)
}

H_gain <- function(P, H, aPH = parlist[["aPH"]]) {
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}

H_loss <- function(H, X, V, M, q, hM = parlist[["hM"]], .no_M_to_X = FALSE) {
    aX <- parlist[["aX"]]
    hX <- parlist[["hX"]]
    if (!.no_M_to_X) {
        Hl <- ((aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / H
    } else {
        Hl <- (aX*H*X)/(1 + aX*hX*(V + H)) / H
    }
    return(Hl)
}





do_sim_fig4 <- function(.midges_not_to) {

    stopifnot(length(.midges_not_to) == 1 && .midges_not_to %in% c("X", "none", "D"))

    q_ <- 3

    sim <- food_web(tmax = 100, s = 10, b = 20, w = 20,
                    other_pars = list(q = 3),
                    midges_not_to = .midges_not_to) %>%
        spread(pool, N) %>%
        mutate(Vg = V_gain(detritivore, detritus),
               Vl = V_loss(detritivore, predator, herbivore, midge, q_,
                           .no_M_to_X = (.midges_not_to == "X")),
               Hg = H_gain(plant, herbivore),
               Hl = H_loss(herbivore, predator, detritivore, midge, q_,
                           .no_M_to_X = (.midges_not_to == "X"))) %>%
        select(-soil:-midge) %>%
        gather("variable", "value", Vg:Hl) %>%
        mutate(pool = factor(ifelse(grepl("^V", variable), "detritivore", "herbivore")),
               type = factor(ifelse(grepl("l$", variable), "top-down", "bottom-up"))) %>%
        select(-variable) %>%
        select(pool, type, everything()) %>%
        group_by(pool, type) %>%
        mutate(value = value - value[1]) %>%
        ungroup() %>%
        mutate(not_to = .midges_not_to)

    return(sim)
}


fig4_df <- tibble(.midges_not_to = c("X", "none", "D")) %>%
    pmap_dfr(do_sim_fig4)

fig4_panel <- function(.midges_not_to, .ylims = c(-0.037, 0.1)) {

    # .midges_not_to = "D"; .ylims = c(-0.037, 0.1)
    # rm(.midges_not_to, .ylims)

    stopifnot(length(.midges_not_to) == 1 && .midges_not_to %in% c("none", "D", "X"))


    .title <- case_when(.midges_not_to == "none" ~ "Midges to detritus and predators",
                        .midges_not_to == "D" ~ "Midges to predators only",
                        .midges_not_to == "X" ~ "Midges to detritus only",
                        TRUE ~ "ERROR")


    .td_pool <- "detritivore"


    mt_combo <- fig4_df %>%
        filter(not_to == .midges_not_to) %>%
        select(-not_to)



    mt_plot <- ggplot(data = NULL) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(data = mt_combo %>% filter(type == "bottom-up"),
                  aes(time, value,
                      color = pool,
                      group = interaction(pool, type)), size = 1, linetype = 1) +
        geom_line(data = mt_combo %>% filter(type == "top-down", pool == .td_pool),
                  aes(time, value), size = 1, color = "gray60") +
        xlab("Time (days)") +
        ggtitle(.title) +
        scale_color_manual(values = color_pal()[1:2]) +
        scale_y_continuous(expression("Effect on pool (" * day^{-1} * ")" ),
                           limits = .ylims, breaks = c(0, 0.05, 0.1)) +
        theme(legend.position = "none",
              strip.text = element_text(face = "plain", size = 10,
                                        margin = margin(b = 4)),
              panel.spacing.x = unit(1.5, "lines"),
              strip.background = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.margin = margin(0,t=12,b=12,r=4),
              plot.title = element_text(size = 12, hjust = 0.5, face = "plain",
                                        margin = margin(0,0,0,b = 10)))


    if (.midges_not_to == "X") {

        mt_plot <- mt_plot +
            geom_text(data = tibble(time =  c(  38,    16,    25),
                                    value = c(0.06, 0.034, -0.02),
                                    lab = sprintf("italic(%s)",
                                                  c("'BU'['detritivore']",
                                                    "'BU'['herbivore']",
                                                    "'TD'['both']"))),
                      aes(time, value, label = lab),
                      color = c(color_pal()[1:2], "gray60"),
                      hjust = 0, vjust = 0, lineheight = 0.75, size = 10 / 2.835,
                      parse = TRUE)

    }

    if (.midges_not_to == "none") {

        mt_plot <- mt_plot +
            theme(axis.title.y = element_text(margin = margin(0,0,0,r=12),
                                              size = 11))

    } else {

        mt_plot <- mt_plot +
            theme(axis.title.y = element_blank())

    }

    if (.midges_not_to != "D") {
        mt_plot <- mt_plot +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    } else {
        mt_plot <- mt_plot +
            geom_segment(data = tibble(x = c(50, 37) + 8,
                                       xend = c(40, 22),
                                       y = c(0.05, -0.03) + c(0, 0.002),
                                       yend = y + c(-0.02, 0.01)),
                         aes(x, y, xend = xend, yend = yend),
                         # curvature = 1,
                         arrow = arrow(length = unit(0.3, "lines"))) +
            geom_text(data = tibble(lab = c("TD[intensification]", "TD[alleviation]"),
                                    time = c(50, 37) + 10,
                                    value = c(0.05, -0.03),
                                    hj = c(0, 0)),
                      aes(time, value, label = lab, hjust = hj),
                      parse = TRUE, vjust = 0.5) +
            theme(axis.title.x = element_text(margin = margin(0,0,0,t=4),
                                              size = 11))
    }

    return(mt_plot)

}



fig4_panel_list <- tibble(.midges_not_to = c("X", "none", "D")) %>%
    pmap(fig4_panel) %>%
    map(~ .x + theme(axis.line.y.left = element_line(size = 1)))



cairo_pdf(file = paste0(dir, "4-up_down_attack_rates.pdf"), width = 3.5, height = 6.5)
ggarrange(plots = fig4_panel_list, ncol = 1, labels = sprintf("(%s)", letters[1:3]),
          label.args = list(gp = gpar(fontsize = 14, fontface =  "plain"),
                            vjust = 2, hjust = -1.75))
dev.off()







# ================================================================================*
# ================================================================================*

# Write captions ----

# ================================================================================*
# ================================================================================*

writeLines(sprintf("Figure %i. %s\n", 2:4, c(fig2_caption, fig3_caption, fig4_caption)),
           paste0(dir, "captions.txt"))
