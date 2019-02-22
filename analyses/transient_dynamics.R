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


middle_sim %>%
    filter(pool != "midge") %>%
    mutate(pool = droplevels(pool),
           level = ifelse(pool %in% c("soil", "plant", "detritus"), 1, 0) %>%
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
    facet_wrap(~ level, nrow = 2, scales = "free_x") +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none") +
    geom_text(data = tibble(
        pool = sort(unique(middle_sim$pool[middle_sim$pool != "midge"])),
        time =  c(43,  17, 80,  60,  85,  50),
        N =     c(2.1, 3.5,  1.8, 1.5, 1.8, 0.8),
        level = factor(rep(1:0, each = 3), levels = 0:1,
                       labels = paste(c("Upper", "Lower"), "trophic levels"))),
        aes(label = pool, color = pool), fontface = "bold") +
    geom_text(data = tibble(time = 22, N = 0.75),
              label = "midge", color = "gray40", fontface = "bold")



