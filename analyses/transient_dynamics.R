# ======================================
# Preliminaries
# ======================================

# load packages
library(mtf)
library(tidyverse)
library(forcats)



middle_sim <- food_web(tmax = 100, s = 10, b = 50, w = 20, .lM = 0.1, .aM = 0.01) %>%
    mutate(pool = fct_recode(pool, soil = "nitrogen"))

# RColorBrewer::brewer.pal(6, "Dark2")

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
    scale_y_continuous("Scaled N") +
    scale_x_continuous("Time") +
    facet_wrap(~ level, nrow = 2) +
    scale_color_brewer(palette = "Dark2") +
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
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_ribbon(data = other_sims %>%
                    filter(pool == "midge") %>%
                    group_by(area, aM) %>%
                    mutate(N = N - N[1]) %>%
                    ungroup() %>%
                    mutate(N = N / sd(N)),
                aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Scaled N") +
    scale_x_continuous("Time") +
    facet_wrap(~ aM + area) +
    scale_color_brewer(NULL, palette = "Dark2") +
    theme(legend.position = "none",
          strip.text.x = element_blank(),
          panel.spacing.y = unit(1.5, "lines")) +
    geom_text(data = tibble(pool = factor(upper_levels),
                            time =  c(40, 30, 40),
                            N =     c(2.6, 1.2, 0.3),
                            aM = factor(rep(1e-4, 3), levels = levels(other_sims$aM)),
                            area = factor(rep(1000, 3), levels = levels(other_sims$area))),
        aes(label = pool, color = pool), fontface = "bold", hjust = 0) +
    geom_text(data = tibble(time =  rep(0, 2),
                            N =     rep(4.44, 2),
                            aM = factor(c(1e-4, 1), levels = levels(other_sims$aM)),
                            area = factor(rep(100, 2), levels = levels(other_sims$area)),
                            lab = c("Low attack rates on midges:",
                                    "High attack rates on midges:")),
              aes(label = lab), hjust = 0, vjust = 1,
              fontface = "bold", size = 12 * (5/14)) +
    coord_cartesian(ylim = c(-1.4121, 4.443)) +
    NULL




