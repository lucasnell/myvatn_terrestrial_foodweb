# ======================================
# Plots for transient dynamics
# ======================================



#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
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



middle_sim <- food_web(tmax = 100, s = 10, b = 50, w = 20,
                       other_pars = list(aM = 0.001)) %>%
    mutate(pool = fct_recode(pool, soil = "nitrogen"))

# < order of colors: green, red, purple, pink, light green, yellow >
# rows correspond to `RColorBrewer::brewer.pal(6, "Dark2")`
rgb_mat <- rbind(c(27,158,119), c(217,95,2), c(117,112,179),
                 c(231,41,138), c(102,166,30), c(230,171,2))
# Switching detritus with plant and detritivore with herbivore
# Order is now "detritivore", "herbivore", "predator", "soil", "detritus", "plant"
rgb_mat <- rgb_mat[c(2, 5, 3, 6, 4, 1),]
# Multipliers for each color, < 1 makes it darker
rgb_mults <- c(0.5, 1.3, 1.1,
               1.05, 0.9, 1)
color_pal <- apply(rgb_mat * matrix(rgb_mults, 6, 3), 1,
                   function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))





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
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_ribbon(data = middle_sim %>%
                    filter(pool == "midge") %>%
                    mutate(N = N / sd(N)),
                aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Proportional change in nitrogen", breaks = c(0, 1, 2)) +
    scale_x_continuous("Time (days)") +
    facet_wrap(~ level, nrow = 2) +
    scale_color_manual(values = color_pal[c(4:6, 1:3)]) +
    theme(legend.position = "none",
          panel.spacing = unit(1.5, "lines")) +
    geom_text(data = tibble(
        pool = sort(unique(middle_sim$pool[middle_sim$pool != "midge"])),
        time =  c( 43,  50,  50,  45,  50,  60),
        N =     c(2.1, 0.8, 0.2, 0.5, 0.1, 1.0),
        level = factor(rep(1:0, each = 3), levels = 0:1,
                       labels = paste(c("Upper", "Lower"), "trophic levels"))),
        aes(label = pool, color = pool), size = 10 / 2.835) +
    geom_text(data = tibble(time = 22, N = 2,
                            level = factor(0, levels = 0:1,
                                           labels = paste(c("Upper", "Lower"),
                                                          "trophic levels"))),
              label = "midge", color = "gray40", size = 10 / 2.835) +
    geom_text(data = tibble(time =  rep(0, 2), N = rep(2.523417, 2),
                            level = factor(paste(c("Upper", "Lower"), "trophic levels")),
                            labs = letters[1:2]),
              aes(label = labs), hjust = 0, vjust = 1, size = 12 / 2.835)
# trans_p1


# ggsave("~/Desktop/2-N_timeseries.pdf", trans_p1, width = 5, height = 5)





other_sims <- map2_dfr(rep(c(100, 1000), each = 2),
                       rep(c(1e-4, 1), 2),
                       function(area_, aM_) {
                           b_ <- area_ / 20
                           food_web(tmax = 100, s = 10, b = b_, w = 20,
                                    other_pars = list(aM = aM_, lM = 0.1)) %>%
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
    mutate(N = (N - N[1]) / N[1]) %>%
    ungroup() %>%
    mutate(aM = factor(aM, levels = range(as.numeric(paste(aM))),
                       labels = paste(c("low", "high"), "attack rate")),
           area = factor(area, levels = range(as.numeric(paste(area))),
                      labels = paste(c("low", "high"), "midge input"))) %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_ribbon(data = other_sims %>%
                    filter(pool == "midge") %>%
                    group_by(area, aM) %>%
                    mutate(N = N - N[1]) %>%
                    ungroup() %>%
                    mutate(N = N / sd(N),
                           aM = factor(aM, levels = range(as.numeric(paste(aM))),
                                       labels = paste(c("low", "high"), "attack rate")),
                           area = factor(area, levels = range(as.numeric(paste(area))),
                                         labels = paste(c("low", "high"), "midge input"))),
                aes(ymin = 0, ymax = N), fill = "gray80", color = NA) +
    geom_line(aes(color = pool), size = 0.75) +
    scale_y_continuous("Proportional change in nitrogen") +
    scale_x_continuous("Time (days)") +
    scale_color_manual(NULL, values = color_pal) +
    geom_text(data = tibble(time =  rep(0, 4), N = rep(4.4, 4),
                            aM = factor(paste(rep(c("low", "high"), each=2), "attack rate"),
                                        levels = paste(c("low", "high"), "attack rate")),
                            area = factor(paste(rep(c("low", "high"), 2), "midge input"),
                                          levels = paste(c("low", "high"), "midge input")),
                            labs = letters[1:4]),
        aes(label = labs), hjust = 0, vjust = 1, size = 12 / 2.835) +
    geom_text(data = tibble(pool = factor(upper_levels),
                            time =  c(50, 90, 100),
                            N =     c(1.45, -0.15, 0.97),
                            aM = factor(paste(rep("low", 3), "attack rate"),
                                        levels = paste(c("low", "high"), "attack rate")),
                            area = factor(paste(rep("high", 3), "midge input"),
                                          levels = paste(c("low", "high"), "midge input"))),
        aes(label = pool, color = pool), hjust = 1, size = 10 / 2.835) +
    facet_grid(aM ~ area) +
    theme(legend.position = "none",
          strip.text.y = element_text(face = "plain", size = 11, angle = 270,
                                      margin = margin(l = 4)),
          strip.text.x = element_text(face = "plain", size = 11,
                                      margin = margin(b = 4)),
          panel.spacing.y = unit(1.5, "lines"),
          panel.spacing.x = unit(2.5, "lines"),
          axis.title = element_text(size = 12)) +
    NULL
# trans_p2

# ggsave("~/Desktop/3-N_midge_attack.pdf", trans_p2, width = 5, height = 5)





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
                                    other_pars = list(aM = aM_, lM = 0.1)) %>%
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
           aM = factor(aM, levels = range(as.numeric(paste(aM))),
                       labels = paste(c("low", "high"), "attack rate")),
           area = factor(area, levels = range(as.numeric(paste(area))),
                         labels = paste(c("low", "high"), "midge input"))) %>%
    select(-variable) %>%
    select(area, aM, pool, type, everything()) %>%
    group_by(area, aM, pool, type) %>%
    mutate(value = value - value[1]) %>%
    ungroup()







trans_p3 <- ggplot(data = NULL) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_line(data = other_sims2 %>% filter(type == "bottom-up"),
              aes(time, value,
                  # linetype = type,
                  color = pool,
                  group = interaction(pool, type)), size = 1, linetype = 1) +
    geom_line(data = other_sims2 %>% filter(type == "top-down", pool == "detritivore"),
                  aes(time, value,
                      # linetype = type,
                      # color = pool,
                      group = interaction(pool, type)), size = 1, color = "gray60") +
    geom_text(data = tibble(time =  rep(0, 4), N = rep(max(other_sims2$value), 4),
                            aM = factor(paste(rep(c("low", "high"), each=2), "attack rate"),
                                        levels = paste(c("low", "high"), "attack rate")),
                            area = factor(paste(rep(c("low", "high"), 2), "midge input"),
                                          levels = paste(c("low", "high"), "midge input")),
                            labs = letters[1:4]),
              aes(time, N, label = labs), hjust = 0, vjust = 1, size = 12 / 2.835) +
    scale_y_continuous(expression("Effect on pool (" * day^{-1} * ")" )) +
    scale_x_continuous("Time (days)") +
    scale_color_manual(values = color_pal[1:2]) +
    geom_text(data = tibble(time =  c(  40,   18,    23),
                            value = c(0.10, 0.055, 0.022),
                            aM = factor(paste(c("low","low","low"), "attack rate"),
                                        levels = levels(other_sims2$aM)),
                            area = factor(paste(c("high","high","high"), "midge input"),
                                          levels = levels(other_sims2$area)),
                            lab = c("BU detritivore", "BU\nherbivore",  "TD both")),
              aes(time, value, label = lab), color = c(color_pal[1:2], "gray60"),
              hjust = 0, vjust = 0, lineheight = 0.75, size = 10 / 2.835) +
    facet_grid(aM ~ area) +
    theme(legend.position = "none",
          strip.text.y = element_text(face = "plain", size = 11, angle = 270,
                                      margin = margin(l = 4)),
          strip.text.x = element_text(face = "plain", size = 11,
                                      margin = margin(b = 4)),
          panel.spacing.y = unit(1.5, "lines"),
          panel.spacing.x = unit(2.5, "lines"),
          axis.title = element_text(size = 12)) +
    NULL




# ggsave("~/Desktop/4-up_down_attack_rates.pdf", trans_p3, width = 5, height = 5)






# trans_p1
# trans_p2
# trans_p3
