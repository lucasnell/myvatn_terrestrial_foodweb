# ======================================
# Plots for transient dynamics
# ======================================



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
})



dir <- sprintf("~/Box Sync/Iceland Food Web Model/Results/Figures_%s/", Sys.Date())

if (!dir.exists(dir)) dir.create(dir)


middle_sim <- food_web(tmax = 100, s = 10, b = 50, w = 20,
                       other_pars = list(f = 3))


upper_levels <- c("detritivore", "herbivore", "predator")


trans_p1 <- middle_sim %>%
    filter(pool != "midge") %>%
    mutate(pool = droplevels(pool),
           level = ifelse(pool %in% upper_levels, 0, 1) %>%
               factor(levels = 0:1,
                      labels = paste(c("Upper", "Lower"), "trophic levels"))) %>%
    group_by(pool) %>%
    mutate(N = (N - N[1]) / (N[1])) %>%
    ungroup() %>%
    {max_N <<- max(.$N); .} %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, color = "gray70") +
    # geom_ribbon(data = middle_sim %>%
    #                 filter(pool == "midge") %>%
    #                 mutate(N = N / sd(N)),
    #             aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_segment(data = tibble(time = 10, time2 = 10+20, N = -0.25),
                 aes(xend = time2, yend = N), size = 1.5) +
    # geom_errorbarh(data = tibble(time = 10, time2 = 10+20, N = 2),
    #                aes(xmin = time, xmax = time2), height = 0.25, size = 1) +
    geom_line(aes(color = pool), size = 1) +
    scale_y_continuous("Proportional change in N", breaks = c(0, 2, 4),
                       limits = c(-0.5, max_N)) +
    scale_x_continuous("Time (days)") +
    facet_wrap(~ level, nrow = 2) +
    scale_color_manual(values = color_pal()[c(4:6, 1:3)]) +
    theme(legend.position = "none",
          panel.spacing = unit(1.5, "lines")) +
    geom_text(data = tibble(
        pool = sort(unique(middle_sim$pool[middle_sim$pool != "midge"])),
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
    geom_text(data = tibble(time =  rep(0, 2),
                            N = rep(max_N, 2),
                            level = factor(paste(c("Upper", "Lower"), "trophic levels")),
                            labs = letters[1:2]),
              aes(label = labs), hjust = 0, vjust = 1, size = 12 / 2.835)
# trans_p1


ggsave(paste0(dir, "2-N_timeseries.pdf"), trans_p1, width = 5, height = 5)





other_sims <- map_dfr(c(8e-3, 8),
                       function(f_) {
                           food_web(tmax = 100, s = 10, b = 20, w = 20,
                                    other_pars = list(f = f_)) %>%
                               filter(pool %in% c(upper_levels, "midge")) %>%
                               mutate(pool = droplevels(pool),
                                      f = f_)
    }) %>%
    mutate_at(vars(f), factor) %>%
    mutate(pool = factor(paste(pool), levels = c(upper_levels, "midge")))


trans_p2 <- other_sims %>%
    filter(!pool %in% c("midge", "predator")) %>%
    mutate(pool = droplevels(pool)) %>%
    group_by(pool, f) %>%
    mutate(N = (N - N[1]) / N[1]) %>%
    ungroup() %>%
    mutate(f = factor(f, levels = range(as.numeric(paste(f))),
                       labels = paste(c("low", "high"), "accessibility"))) %>%
    {max_N <<- max(.$N) * (2 / 1.5); .} %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_segment(data = tibble(time = 10, time2 = 10+20, N = -0.05),
                 aes(xend = time2, yend = N), size = 1.5) +
    geom_ribbon(data = other_sims %>%
                    filter(pool == "predator") %>%
                    group_by(f) %>%
                    mutate(N = N - N[1]) %>%
                    ungroup() %>%
                    mutate(N = N / 1.5,
                           f = factor(f, levels = range(as.numeric(paste(f))),
                                       labels = paste(c("low", "high"), "accessibility"))),
                aes(ymin = 0, ymax = N), fill = color_pal(0.5)[3]) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Proportional change in N", limits = c(-0.15, max_N),
                       sec.axis = sec_axis(~ . * 1.5,
                                           name = "Proportional change in N\n(predator)",
                                           breaks = 0:2)) +
    scale_x_continuous("Time (days)") +
    scale_color_manual(NULL, values = color_pal()) +
    geom_text(data = tibble(time =  rep(0, 2), N = rep(max_N, 2),
                            f = factor(paste(c("low", "high"), "accessibility"),
                                        levels = paste(c("low", "high"), "accessibility")),
                            labs = letters[1:2]),
        aes(label = labs), hjust = 0, vjust = 1, size = 12 / 2.835) +
    geom_text(data = tibble(time =  75,
                            N =     1.4,
                            f = factor("high accessibility",
                                        levels = paste(c("low", "high"), "accessibility"))),
        label = "predator", color = color_pal()[3], hjust = 1, vjust = 1,
        size = 10 / 2.835) +
    geom_text(data = tibble(pool = factor(upper_levels[upper_levels != "predator"]),
                            time =  c(92, 100),
                            N =     c(0.3, -0.1),
                            f = factor(paste(rep("high", 2), "accessibility"),
                                        levels = paste(c("low", "high"), "accessibility"))),
        aes(label = pool, color = pool), hjust = 1, size = 10 / 2.835) +
    geom_text(data = tibble(time = 20, N = -0.1,
                            f = factor("high accessibility",
                                       levels = paste(c("low", "high"), "accessibility"))),
              label = "pulse", size = 10 / 2.835, hjust = 0.5, vjust = 1,
              color = "black") +
    facet_grid( ~ f) +
    theme(legend.position = "none",
          strip.text.y = element_text(face = "plain", size = 11, angle = 270,
                                      margin = margin(l = 4)),
          strip.text.x = element_text(face = "plain", size = 11,
                                      margin = margin(b = 4)),
          panel.spacing.y = unit(1.5, "lines"),
          panel.spacing.x = unit(2.5, "lines")) +
    NULL
# trans_p2

ggsave(paste0(dir, "3-N_midge_attack.pdf"), trans_p2, width = 5, height = 3)




parlist <- par_estimates %>%
    filter(V==1, H==1, R==1, iN == 10) %>%
    as.list()
V_gain <- function(V, D) {
    aDV <- parlist[["aDV"]]
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
V_loss <- function(V, R, H, M, f) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    hM <- parlist[["hM"]]
    ((aR*V*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / V
}
H_gain <- function(P, H) {
    aPH <- parlist[["aPH"]]
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}
H_loss <- function(H, R, V, M, f) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    hM <- parlist[["hM"]]
    ((aR*H*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / H
}




other_sims2 <- map_dfr(c(8e-3, 8),
                       function(f_) {
                           food_web(tmax = 100, s = 10, b = 20, w = 20,
                                    other_pars = list(f = f_)) %>%
                               mutate(f = f_)
                       }) %>%
    spread(pool, N) %>%
    mutate(Vg = V_gain(detritivore, detritus),
           Vl = V_loss(detritivore, predator, herbivore, midge, f),
           Hg = H_gain(plant, herbivore),
           Hl = H_loss(herbivore, predator, detritivore, midge, f)) %>%
    select(-soil:-midge) %>%
    gather("variable", "value", Vg:Hl) %>%
    mutate(pool = factor(ifelse(grepl("^V", variable), "detritivore", "herbivore")),
           type = factor(ifelse(grepl("l$", variable), "top-down", "bottom-up")),
           f = factor(f, levels = range(as.numeric(paste(f))),
                       labels = paste(c("low", "high"), "accessibility"))) %>%
    select(-variable) %>%
    select(f, pool, type, everything()) %>%
    group_by(f, pool, type) %>%
    mutate(value = value - value[1]) %>%
    ungroup()







trans_p3 <- ggplot(data = NULL) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(data = other_sims2 %>% filter(type == "bottom-up"),
              aes(time, value,
                  color = pool,
                  group = interaction(pool, type)), size = 1, linetype = 1) +
    geom_line(data = other_sims2 %>% filter(type == "top-down", pool == "detritivore"),
                  aes(time, value), size = 1, color = "gray60") +
    geom_text(data = tibble(time =  rep(0, 2), N = rep(max(other_sims2$value), 2),
                            f = factor(paste(c("low", "high"), "accessibility"),
                                        levels = paste(c("low", "high"), "accessibility")),
                            labs = letters[1:2]),
              aes(time, N, label = labs), hjust = 0, vjust = 1, size = 12 / 2.835) +
    scale_y_continuous(expression("Effect on pool (" * day^{-1} * ")" ),
                       limits = c(-0.03615922, 0.08581619)) +
    scale_x_continuous("Time (days)") +
    scale_color_manual(values = color_pal()[1:2]) +
    geom_text(data = tibble(time =  c(  40,   17,    35),
                            value = c(0.07, 0.033, 0.003),
                            f = factor(paste(c("low","low","low"), "accessibility"),
                                        levels = levels(other_sims2$f)),
                            lab = c("BU detritivore", "BU\nherbivore",  "TD both")),
              aes(time, value, label = lab), color = c(color_pal()[1:2], "gray60"),
              hjust = 0, vjust = 0, lineheight = 0.75, size = 10 / 2.835) +
    facet_grid(~ f) +
    theme(legend.position = "none",
          strip.text = element_text(face = "plain", size = 11,
                                    margin = margin(b = 4)),
          panel.spacing.x = unit(2.5, "lines"),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0, size = 12),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12),
          legend.margin = margin(0,0,0,0)) +
    NULL
# trans_p3




impute <- function(time, value) {
    X_ <- zoo::zoo(value, time) %>%
        zoo::na.approx() %>%
        as.numeric()
    return(X_)
}




trans_p4 <- other_sims2 %>%
    filter(type == "top-down", pool == "detritivore", f == "high accessibility") %>%
    select(-pool, -f, -type) %>%
    {bind_rows(., tibble(time = seq(34.1, 34.9, 0.1), value = NA_real_))} %>%
    arrange(time) %>%
    mutate(value = impute(time, value)) %>%
    {
        ggplot(., aes(time, value)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_ribbon(data = . %>% filter(time < 35),
                        aes(ymin = 0, ymax = value),
                        fill = "gray80") +
            geom_ribbon(data = . %>% filter(time >= 35),
                        aes(ymin = 0, ymax = value),
                        fill = "gray80") +
            geom_line(size = 1, color = "gray60") +
            geom_text(data = tibble(time = c(36, 14),
                             value = c(0.015, -0.02),
                             lab = paste0("top-down\n", c("intensification",
                                                          "alleviation"))),
                      aes(label = lab), hjust = 0, vjust = 0.5, size = 10 / 2.835) +
            geom_text(data = tibble(time =  0, N = max(other_sims2$value),
                                    labs = letters[3]),
                      aes(time, N, label = labs), hjust = 0, vjust = 1,
                      size = 12 / 2.835) +
            scale_y_continuous(expression("Effect on pool (" * day^{-1} * ")" ),
                               limits = c(-0.03615922, 0.08581619)) +
            scale_x_continuous("Time (days)") +
            theme(plot.title = element_text(hjust = 0, size = 12),
                  axis.text = element_text(size = 10, color = "black"),
                  axis.title = element_text(size = 12),
                  legend.margin = margin(0,0,0,0),
                  legend.text = element_text(size = 10),
                  strip.background = element_blank())
    }




pdf(file = paste0(dir, "4-up_down_attack_rates.pdf"), width = 5, height = 6)
cowplot::plot_grid(trans_p3, trans_p4, ncol = 1)
dev.off()


# trans_p1
# trans_p2
# trans_p3
# trans_p4
