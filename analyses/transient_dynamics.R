# ======================================
# Plots for transient dynamics
# ======================================


#'
#' ----------
#' To do list:
#' ----------
#'
#' * Sort out figure sizes
#' * Adjust labels for colors once sizes are decided
#' * Sizes of axis labels
#'

# trans_p1
# trans_p2
# trans_p3



# load packages
library(mtf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(forcats)



middle_sim <- food_web(tmax = 100, s = 10, b = 50, w = 20, .lM = 0.1, .aM = 0.01) %>%
    mutate(pool = fct_recode(pool, soil = "nitrogen"))

# < order of colors: pink, light green, yellow, green, red, purple >
color_pal <- RColorBrewer::brewer.pal(6, "Dark2")
# Switching detritus with plant and detritivore with herbivore
color_pal <- color_pal[c(2, 1, 3, 4, 6, 5)]


upper_levels <- c("detritivore", "herbivore", "predator")


trans_p1 <- middle_sim %>%
    filter(pool != "midge") %>%
    mutate(pool = droplevels(pool),
           level = ifelse(pool %in% upper_levels, 0, 1) %>%
               factor(levels = 0:1, labels = paste(c("Upper", "Lower"), "trophic levels"))) %>%
    group_by(pool) %>%
    mutate(N = (N - N[1]) / sd(N)) %>%
    ungroup() %>%
    ggplot(aes(time, N)) +
    geom_ribbon(data = middle_sim %>%
                    filter(pool == "midge") %>%
                    mutate(N = (N - N[1]) / sd(N)),
                aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Scaled nitrogen") +
    scale_x_continuous("Time (days)") +
    facet_wrap(~ level, nrow = 2) +
    scale_color_manual(values = color_pal[c(4:6, 1:3)]) +
    theme(legend.position = "none",
          strip.text = element_text(size = 12, margin = margin(b = 1, t = 0, 0, 0)),
          panel.spacing = unit(1.5, "lines")) +
    geom_text(data = tibble(
        pool = sort(unique(middle_sim$pool[middle_sim$pool != "midge"])),
        time =  c(43,  17, 80,  60,  85,  50),
        N =     c(2.1, 3.5,  1.8, 1.5, 1.8, 0.8),
        level = factor(rep(1:0, each = 3), levels = 0:1,
                       labels = paste(c("Upper", "Lower"), "trophic levels"))),
        aes(label = pool, color = pool), fontface = "bold") +
    geom_text(data = tibble(time = 22, N = 0.75),
              label = "midge", color = "gray40", fontface = "bold")



# ggsave("~/Desktop/N_timeseries.pdf", trans_p1, width = 5, height = 5)





other_sims <- map2_dfr(rep(c(100, 1000), each = 2), rep(c(1e-4, 1), 2),
                       function(area_, aM_) {
                           b_ <- area_ / 20
                           food_web(tmax = 100, s = 10, b = b_, w = 20,
                                    .lM = 0.1, .aM = aM_) %>%
                               filter(pool %in% c(upper_levels, "midge")) %>%
                               mutate(pool = droplevels(pool),
                                      area = area_, aM = aM_)
    }) %>%
    mutate_at(vars(area, aM), factor) %>%
    mutate(pool = factor(paste(pool), levels = c(upper_levels, "midge")))


trans_p2 <- other_sims %>%
    filter(pool != "midge") %>%
    mutate(pool = droplevels(pool)) %>%
    group_by(pool, area, aM) %>%
    mutate(N = N - N[1]) %>%
    group_by(pool) %>%
    mutate(N = N / sd(N)) %>%
    ungroup() %>%
    mutate(aM = factor(aM, levels = c(1e-4, 1), labels = c("low", "high"))) %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_ribbon(data = other_sims %>%
                    filter(pool == "midge") %>%
                    group_by(area, aM) %>%
                    mutate(N = N - N[1]) %>%
                    ungroup() %>%
                    mutate(N = N / sd(N),
                           aM = factor(aM, levels = c(1e-4, 1), labels = c("low", "high"))),
                aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Scaled nitrogen") +
    scale_x_continuous("Time (days)") +
    scale_color_manual(NULL, values = color_pal) +
    geom_text(data = tibble(pool = factor(upper_levels),
                            time =  c(40, 26.5, 40),
                            N =     c(2.6, 1.2, 0.3),
                            aM = factor(rep("low", 3), levels = c("low", "high")),
                            area = factor(rep(1e3, 3), levels = levels(other_sims$area))),
        aes(label = pool, color = pool), fontface = "bold", hjust = 0, size = 3) +
    coord_cartesian(ylim = c(-1.4121, 4.443)) +
    facet_grid(area ~ aM) +
    labs(subtitle = "Attack rate on midges") +
    theme(legend.position = "none",
          strip.text.y = element_blank(),
          panel.spacing.y = unit(1.5, "lines"),
          panel.spacing.x = unit(2.5, "lines"),
          plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5)) +
    NULL


ggsave("~/Desktop/N_midge_attack.pdf", trans_p2, width = 5, height = 5)





V_gain <- function(V, D) {
    parlist <- par_estimates %>%
        filter(V==1, H==1, R==1, iN == formals(food_web)$.iN) %>%
        as.list()
    aDV <- parlist[["aDV"]]
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
V_loss <- function(V, R, H, M, aM) {
    parlist <- par_estimates %>%
        filter(V==1, H==1, R==1, iN == formals(food_web)$.iN) %>%
        as.list()
    aR <- parlist[["aR"]]
    hVHM <- parlist[["hVHM"]]
    ((aR*V*R)/(1 + aR*hVHM*(V + H) + aM*hVHM*M)) / V
}
H_gain <- function(P, H) {
    parlist <- par_estimates %>%
        filter(V==1, H==1, R==1, iN == formals(food_web)$.iN) %>%
        as.list()
    aPH <- parlist[["aPH"]]
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}
H_loss <- function(H, R, V, M, aM) {
    parlist <- par_estimates %>%
        filter(V==1, H==1, R==1, iN == formals(food_web)$.iN) %>%
        as.list()
    aR <- parlist[["aR"]]
    hVHM <- parlist[["hVHM"]]
    ((aR*H*R)/(1 + aR*hVHM*(V + H) + aM*hVHM*M)) / H
}




other_sims2 <- map2_dfr(rep(c(100, 1000), each = 2), rep(c(1e-4, 1), 2),
                       function(area_, aM_) {
                           b_ <- area_ / 20
                           food_web(tmax = 100, s = 10, b = b_, w = 20,
                                    .lM = 0.1, .aM = aM_) %>%
                               mutate(area = area_, aM = aM_)
                       }) %>%
    spread(pool, N) %>%
    mutate(Vg = V_gain(detritivore, detritus),
           Vl = V_loss(detritivore, predator, herbivore, midge, aM),
           Hg = H_gain(plant, herbivore),
           Hl = H_loss(herbivore, predator, detritivore, midge, aM)) %>%
    select(-nitrogen:-midge) %>%
    mutate_at(vars(area), factor) %>%
    gather("variable", "value", Vg:Hl) %>%
    mutate(pool = factor(ifelse(grepl("^V", variable), "detritivore", "herbivore")),
           type = factor(ifelse(grepl("l$", variable), "top-down", "bottom-up")),
           aM = factor(aM, levels = c(1e-4, 1), labels = c("low", "high"))) %>%
    select(-variable) %>%
    select(area, aM, pool, type, everything()) %>%
    group_by(area, aM, pool, type) %>%
    mutate(value = value - value[1])





trans_p3 <- other_sims2 %>%
    ungroup() %>%
    mutate(lty = ifelse(pool == "herbivore" & type == "top-down", 0, 1) %>%
               factor()) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_ribbon(data = other_sims %>%
                    filter(pool == "midge") %>%
                    group_by(area, aM) %>%
                    mutate(N = N - N[1]) %>%
                    ungroup() %>%
                    mutate(N = 0.15 * N / max(N),
                           aM = factor(aM, levels = c(1e-4, 1), labels = c("low", "high"))),
                aes(time, ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(time, value, linetype = lty, color = pool, group = interaction(pool, type)), size = 1) +
    scale_y_continuous(expression("Effect on detritivores (" * day^{-1} * ")" )) +
    scale_x_continuous("Time (days)") +
    facet_grid(area ~ aM) +
    labs(subtitle = "Attack rate on midges") +
    theme(legend.position = "none",
          strip.text.y = element_blank(),
          panel.spacing.y = unit(1.5, "lines"),
          panel.spacing.x = unit(2.5, "lines"),
          plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5)) +
    scale_color_manual(values = color_pal[1:2]) +
    geom_text(data = tibble(time =  c(50, 50),
                            value = c(0.098, 0.00),
                            aM = factor(rep("high", 2), levels = levels(other_sims2$aM)),
                            area = factor(rep(1e3, 2), levels = levels(other_sims2$area)),
                            lab = c("bottom-up", "top-down")),
              aes(time, value, label = lab), fontface = "bold", hjust = 0) +
    scale_linetype_manual(values = c(2, 1)) +
    NULL

# ADD RIGHT-HAND LABEL FOR MIDGE-PULSE INTENSITY
# DO THE SAME THING WITH trans_p2


ggsave("~/Desktop/up_down_attack_rates.pdf", trans_p3, width = 5, height = 5)






# trans_p1
# trans_p2
# trans_p3
