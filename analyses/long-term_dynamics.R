# ======================================
# Preliminaries
# ======================================

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(forcats)
})

pulse_df <- read_csv("data-raw/pulse_data.csv",
                     col_types = cols(
                         w = "f",
                         b = "f",
                         aM = "f",
                         pool = "f",
                         .default = "d"))




heat_plot <- function(.col, .pool, .df = NULL, .w = NULL, .title = NULL,
                      flip_sign = FALSE, flip_color = FALSE,
                      keep_y = TRUE) {
    if (is.null(.df)) .df <- pulse_df
    if (is.null(.w)) .w <- .df$w %>% paste() %>% as.numeric() %>% median()
    .col <- enquo(.col)
    if (is.null(.title)) .title <- .pool
    .df <- .df %>%
        filter(pool == .pool, w == .w) %>%
        mutate(aM = aM %>% paste() %>% as.numeric())
    if (flip_sign) {
        .df <- .df %>% mutate_at(vars(!!.col), ~ . * -1)
    }
    midpt <- .df %>% select(!!.col) %>% .[[1]] %>% range(na.rm = TRUE) %>% median()
    p <- .df %>%
        ggplot(aes(aM, area)) +
        geom_tile(aes(fill = !!.col)) +
        scale_fill_viridis_c(NULL, option = "C", direction = ifelse(flip_color, -1, 1)) +
        ggtitle(.title) +
        ylab("Cumulative midge input") +
        scale_x_continuous("Attack rate on midges", breaks = c(0.005, 0.01)) +
        theme(legend.position = "top", legend.text = element_blank(),
              legend.title = element_text(size = 10),
              legend.key.width = grid::unit(0.02, "npc"),
              plot.title = element_text(hjust = 0.5))
    if (!keep_y) p <- p + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank())
    return(p)
}





#'
#' - Overall label for rows in first two plots?
#' - Labeling herbivores in herbivore-only plot?
#'
#'

hm1 <- heat_plot(cum_pos_loss_V, "detritivore", .title = "Top-down elevation")
hm2 <- heat_plot(cum_neg_loss_V, "detritivore",
                 .title = "Top-down depression", flip_color = TRUE, keep_y = FALSE)
hm3 <- heat_plot(cum_gain_V, "detritivore", .title = "Bottom-up",
                 flip_color = TRUE, keep_y = FALSE)
# heat_plot(cum_gain_H, "herbivore", .title = "Enhancement of bottom-up effect", flip_color = TRUE)


hm <- cbind(ggplotGrob(hm1), ggplotGrob(hm2), ggplotGrob(hm3), size = "last")

grid.newpage()
grid.draw(hm)


pdf(file = "~/Desktop/heatmaps.pdf", width = 7, height = 3.5)
grid.newpage()
grid.draw(hm)
dev.off()


heat_plot(max_loss_H, "herbivore")
heat_plot(min_loss_H, "herbivore")
heat_plot(max_gain_H, "herbivore")




#' Gains are never below starting values
# heat_plot(min_gain_V, "detritivore")
# heat_plot(min_gain_H, "herbivore")


pulse_df


p1 <- pulse_df %>%
    filter(pool == "detritivore",
           w %in% c(10,30)) %>%
    mutate_at(vars(aM, b), function(x) x %>% paste() %>% as.numeric()) %>%
    # ggplot(aes(max_loss_V, min_loss_V)) +
    ggplot(aes(cum_pos_loss_V, cum_neg_loss_V)) +
    geom_abline(slope = -1, intercept = 0, linetype = 2, color = "gray50") +
    # geom_point(aes(color = aM)) +
    geom_path(aes(group = interaction(b, factor(w))), color = "gray40") +
    geom_point(aes(color = b), size = 1) +
    # geom_point(color = "dodgerblue", alpha = 0.5) +
    # scale_color_gradient2("Attack rate",
    # scale_color_gradient2("Cum. input",
    #                      low = "dodgerblue", mid = "white", high = "firebrick",
    #                      # midpoint = 0.0055) +
    #                      midpoint = 1500) +
    scale_color_viridis_c("Midge intensity", option = "D", begin = 0.0, end = 1.0) +
    # annotate("text", x = 0.04, y = -0.005, label = "30", hjust = 0, vjust = 0.5) +
    # annotate("text", x = 0.018, y = -0.0075, label = "10", hjust = 1, vjust = 0.5) +
    # ggtitle("detritivore") +
    # xlab("Enhancement of top-down effect") +
    # ylab("Alleviation of top-down effect") +
    xlab("Total elevated top-down effect") +
    ylab("Total depressed top-down effect") +
    # facet_wrap(~ w) +
    theme(legend.position = "top", legend.text = element_blank(),
          legend.title = element_text(size = 10),
          legend.key.width = grid::unit(0.02, "npc"))

p2 <- pulse_df %>%
    filter(pool == "detritivore",
           w %in% c(10,30)) %>%
    mutate_at(vars(aM, b), function(x) x %>% paste() %>% as.numeric()) %>%
    # ggplot(aes(max_loss_V, min_loss_V)) +
    ggplot(aes(cum_pos_loss_V, cum_neg_loss_V)) +
    geom_abline(slope = -1, intercept = 0, linetype = 2, color = "gray50") +
    geom_path(aes(group = interaction(b, factor(w))), color = "gray40") +
    geom_point(aes(color = aM), size = 1) +
    scale_color_viridis_c("Attack rate", option = "B", begin = 0.2, end = 0.9) +
    # xlab("Enhancement of top-down effect") +
    # ylab("Alleviation of top-down effect") +
    xlab("Total elevated top-down effect") +
    ylab("Total depressed top-down effect") +
    theme(legend.position = "top", legend.text = element_blank(),
          legend.title = element_text(size = 10),
          legend.key.width = grid::unit(0.02, "npc"))



library(grid)

p12 <- cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")

grid.newpage()
grid.draw(p12)


pdf(file = "~/Desktop/elev_depr_topdown.pdf", width = 6, height = 4)
grid.newpage()
grid.draw(p12)
dev.off()




p3 <- pulse_df %>%
    filter(w == 20, pool %in% c("detritivore", "herbivore")) %>%
    group_by(aM, b, pool) %>%
    summarize(cum_pos_loss = ifelse(pool[1] == "detritivore", cum_pos_loss_V[1],
                                    cum_pos_loss_H[1]),
              cum_neg_loss = ifelse(pool[1] == "detritivore", cum_neg_loss_V[1],
                                    cum_neg_loss_H[1]),
              cum_loss = cum_pos_loss + cum_neg_loss,
              cum_gain = ifelse(pool[1] == "detritivore", cum_gain_V[1],
                                cum_gain_H[1])) %>%
    ungroup() %>%
    mutate_at(vars(aM, b), function(x) x %>% paste() %>% as.numeric()) %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_point(aes(color = b), size = 1) +
    scale_color_viridis_c("Midge intensity", option = "D", begin = 0.0, end = 1.0) +
    xlab("Total bottom-up effect") +
    ylab("Total top-down effect") +
    # facet_wrap(~ w) +
    facet_grid(pool ~ .) +
    theme(legend.position = "top", legend.text = element_blank(),
          legend.title = element_text(size = 10),
          legend.key.width = grid::unit(0.02, "npc"),
          strip.text = element_blank()) +
    NULL

p4 <- pulse_df %>%
    filter(w == 20, pool %in% c("detritivore", "herbivore")) %>%
    group_by(aM, b, pool) %>%
    summarize(cum_pos_loss = ifelse(pool[1] == "detritivore", cum_pos_loss_V[1],
                                    cum_pos_loss_H[1]),
              cum_neg_loss = ifelse(pool[1] == "detritivore", cum_neg_loss_V[1],
                                    cum_neg_loss_H[1]),
              cum_loss = cum_pos_loss + cum_neg_loss,
              cum_gain = ifelse(pool[1] == "detritivore", cum_gain_V[1],
                                cum_gain_H[1])) %>%
    ungroup() %>%
    mutate_at(vars(aM, b), function(x) x %>% paste() %>% as.numeric()) %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_point(aes(color = aM), size = 1) +
    scale_color_viridis_c("Attack rate", option = "B", begin = 0.2, end = 0.9) +
    xlab("Total bottom-up effect") +
    ylab("Total top-down effect") +
    facet_grid(pool ~ .) +
    theme(legend.position = "top", legend.text = element_blank(),
          legend.title = element_text(size = 10),
          legend.key.width = grid::unit(0.02, "npc"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    NULL


p34 <- cbind(ggplotGrob(p3), ggplotGrob(p4), size = "last")

grid.newpage()
grid.draw(p34)


pdf(file = "~/Desktop/top_vs_bottom.pdf", width = 7, height = 4)
grid.newpage()
grid.draw(p34)
dev.off()





#'
#' Plot info:
#' - return time ~ area
#' - separated by pool and midge leakage
#' - colored by width
#' - lower trophic levels only
#'
pulse_df %>%
    filter(pool %in% c("plant", "soil", "detritus"),
           # mM == 0.5, aM == 5e-3, lM != 0.325,
           aM == levels(aM)[length(levels(aM)) %/% 2],
           # w %in% c(10, 30),
           # area <= 1000) %>%
           TRUE) %>%
    mutate(w = as.numeric(paste(w)),
           pool = factor(paste(pool), levels = c("soil", "detritus", "plant")),
           # lM = fct_recode(lM, low = "0.1", mid = "0.325", high = "0.55") %>%
           # lM = fct_recode(lM, low = "0.1", high = "0.55") %>%
           #     fct_relabel(.fun = function(x) paste0(x, "")),
           # w = factor(w)) %>%
           NULL) %>%
    ggplot(aes(area, return_time)) +
    # geom_hline(yintercept = 75, linetype = 2, color = "gray70") +
    geom_point(aes(color = w), shape = 1, alpha = 0.5, na.rm = TRUE) +
    # geom_line(aes(color = w, linetype = lM), size = 0.75) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20) +
    # facet_grid(pool ~ lM) +
    # scale_color_manual(values = c("firebrick", "dodgerblue")) +
    facet_grid( ~ pool) +
    ylab("Return time (days)") +
    xlab(expression("Cumulative midge input (g" ~ m^2 * ")")) +
    # labs(subtitle = "Midge loss from system") +
    theme(plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5))



#'
#' Plot info:
#' - return time ~ area
#' - separated by pool and midge leakage
#' - colored by width
#' - upper trophic levels only
#'
pulse_df %>%
    filter(pool %in% c("detritivore", "herbivore", "predator"),
           # mM == 0.5, aM == 5e-3) %>%
           aM == levels(aM)[length(levels(aM)) %/% 2]) %>%
    mutate(w = as.numeric(paste(w)), pool = factor(paste(pool)),
           # lM = fct_recode(lM, low = "0.1", mid = "0.325", high = "0.55") %>%
           #     fct_relabel(.fun = function(x) paste(x, ""))) %>%
           NULL) %>%
    ggplot(aes(area, return_time)) +
    # geom_hline(yintercept = 75, linetype = 2, color = "gray70") +
    geom_point(aes(color = w), shape = 1, alpha = 0.5, na.rm = TRUE) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20) +
    # facet_grid(pool ~ aM) +
    facet_grid(~ pool) +
    ylab("Return time (days)") +
    xlab(expression("Cumulative midge input (g" ~ m^2 * ")")) +
    # labs(subtitle = "Midge loss from system") +
    theme(plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5))


#'
#' Plot info:
#' - return time ~ area
#' - separated by pool and attack rate
#' - colored by width
#' - herbivores only
#'
pulse_df %>%
    filter(pool == "herbivore") %>%
    mutate(w = as.numeric(paste(w)), pool = factor(paste(pool))) %>%
    # mutate(w = as.numeric(paste(w)), pool = factor(paste(pool)),
    #        aM = fct_recode(aM, low = "0.001", mid = "0.005", high = "0.01") %>%
    #            fct_relabel(.fun = function(x) paste(x, ""))) %>%
    ggplot(aes(area, min_below)) +
    geom_vline(xintercept = c(500, 850, 1200), linetype = 2, color = "gray70") +
    geom_point(aes(color = w), shape = 1, alpha = 0.5, na.rm = TRUE) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20) +
    facet_grid( ~ aM) +
    ylab("Return time (days)") +
    xlab(expression("Cumulative midge input (g" ~ m^2 * ")")) +
    labs(subtitle = "Attack rate on midges") +
    theme(plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5))



#'
#' Plot info:
#' Trying to show time series for right before, during, and after the discontinuity
#'

timeseries <- expand.grid(area = c(500, 850, 1200),
                          aM = c(1e-3, 5e-3, 1e-2),
                          w = c(10, 30)) %>%
    split(row(.)[,1]) %>%
    map_dfr(function(.row) {
        .area <- .row$area
        .aM <- .row$aM
        .w <- .row$w
        .b <- .area / .w
        fw <- food_web(tmax = 250, s = 10, b = .b, w = .w, .aM = .aM,
                       .mM = 0.5, .lM = 0.325)
        fw <- fw %>%
            filter(pool != "midge") %>%
            mutate(pool = droplevels(pool),
                   time = time - (10 + .w)) %>%
            mutate(w = .w, aM = .aM, area = .area) %>%
            select(w, aM, area, everything())
        return(fw)
    })


timeseries %>%
    filter(pool == "herbivore") %>%
    mutate_at(vars(w, aM, area), factor) %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = c(0.1, -0.1) + 1, linetype = 3, color = "gray60") +
    geom_hline(yintercept = 1, linetype = 2, color = "gray70") +
    geom_line(aes(color = w), size = 0.75) +
    facet_grid(area ~ aM, label = label_both) +
    scale_color_manual(values = c("firebrick", "dodgerblue")) +
    xlab("Time after pulse ends (days)") +
    ylab("Nitrogen content (g)")








# # ------------------------
# # Plotting functions
# # ------------------------
#
# # Construct a heatmap plot
# heat_plot <- function(.df, .col, .pool) {
#     .col <- enquo(.col)
#     .df <- .df %>%
#         filter(pool == .pool) %>%
#         mutate_at(vars(lM, mM), factor)
#     midpt <- .df %>% select(!!.col) %>% .[[1]] %>% range(na.rm = TRUE) %>% median()
#     color_lab_ <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
#         trimws() %>% tools::toTitleCase()
#     .df %>%
#         ggplot(aes(w, b)) +
#         geom_tile(aes(fill = !!.col)) +
#         geom_line(data = tibble(w = rep(seq(10, 30, length.out = 100), 5),
#                                 a = rep(seq(250, 1250, 250), each = 100),
#                                 b = a / w),
#                   aes(group = factor(a)), size = 1) +
#         scale_fill_gradient2(color_lab_,
#                              low = "dodgerblue", mid = "white", high = "firebrick",
#                              midpoint = midpt) +
#         facet_grid(mM ~ lM, labeller = label_both) +
#         ggtitle(.pool) +
#         coord_cartesian(ylim = c(0, 100), xlim = c(10, 30)) +
#         ylab("Rate of input") +
#         xlab("Duration of input")
# }
# # Plot of area vs another variable
# area_plot <- function(.df, .col, .pools, .ylab = NULL) {
#     .col <- enquo(.col)
#     .df <- .df %>%
#         filter(pool %in% .pools, mM == 0.5, lM == 0.5) %>%
#         # filter(pool %in% .pools, mM == 0.5, aM == 1e-2) %>%
#         mutate(w = as.numeric(paste(w)), pool = factor(paste(pool), levels = .pools))
#     if (is.null(.ylab)) {
#         .ylab <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
#             trimws() %>% tools::toTitleCase()
#     }
#     .df %>%
#         ggplot(aes(area, return_time)) +
#         geom_hline(yintercept = 75, linetype = 2, color = "gray70") +
#         geom_point(aes(color = w), shape = 1, alpha = 0.5, na.rm = TRUE) +
#         scale_color_gradient2("Width",
#                               low = "firebrick", mid = "gray80", high = "dodgerblue",
#                               midpoint = 20) +
#         # ggtitle(.pool) +
#         facet_grid(pool ~ aM) +
#         # facet_grid(pool ~ lM) +
#         ylab(.ylab) +
#         xlab("Area under curve")
# }
#
#
#
# # ------------------------
# # Plots
# # ------------------------
#
#
#
# area_plot(pulse_df, return_time, c("detritivore", "herbivore", "predator"))
# area_plot(pulse_df, return_time, c("plant", "soil", "detritus"))
#
#
#
#
#
# heat_plot(pulse_df, to_max, "detritivore")
# heat_plot(pulse_df, to_min, "detritivore")
# heat_plot(pulse_df, cum_N, "detritivore")
#
#
#
# pulse_df %>%
#     filter(pool == "detritivore", lM == 0.1) %>%
#     # mutate(lM = factor(lM)) %>%
#     # ggplot(aes(w, return_time)) +
#     ggplot(aes(b, return_time)) +
#     # geom_point(aes(color = w)) +
#     # geom_line(aes(color = b, group = factor(b)), size = 1) +
#     geom_line(aes(color = w, group = factor(w)), size = 1) +
#     # facet_grid(~ lM) +
#     # facet_wrap(~ w, nrow = 5) +
#     scale_color_gradient2("Width",
#                           low = "firebrick", mid = "gray80", high = "dodgerblue",
#                           midpoint = 20) +
#     # scale_color_gradient2("Intensity",
#     #                       low = "firebrick", mid = "gray80", high = "dodgerblue",
#     #                       midpoint = 50) +
#     NULL
#
#
# pulse_df %>%
#     filter(pool == "detritivore") %>%
#     mutate_at(vars(lM, mM), factor) %>%
#     ggplot(aes(cum_N, max_above)) +
#     geom_line(aes(color = w, group = factor(b)), size = 1) +
#     facet_grid(mM ~ lM, labeller = label_both) +
#     scale_color_gradient2("Width",
#                           low = "firebrick", mid = "gray80", high = "dodgerblue",
#                           midpoint = 20)
#
#
#
# heat_plot(pulse_df, min_below, "detritivore")
# heat_plot(pulse_df, max_above, "detritivore")
#
# heat_plot(pulse_df, return_time, "detritivore")
#
#
# heat_plot(pulse_df, return_time, "herbivore")
# heat_plot(pulse_df, return_time, "predator")
#
#
#
#
#
#
# area_plot(pulse_df, return_time, "detritivore") +
#     geom_vline(xintercept = 350, linetype = 2)
# area_plot(pulse_df, return_time, "herbivore")  # +
#     # geom_vline(data = tibble(xint = c(200, 900, 1100), lM = factor(c(0.1, 0.5, 0.9))),
#     #           aes(xintercept = xint), linetype = 2)
# area_plot(pulse_df, return_time, "predator")
#
# area_plot(pulse_df, return_time, "nitrogen")
# area_plot(pulse_df, return_time, "detritus")
# area_plot(pulse_df, return_time, "plant")

